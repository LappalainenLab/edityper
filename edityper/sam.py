#!/usr/bin/env python

"""Create a SAM file from alignment data for EdiTyper"""

from __future__ import division
from __future__ import print_function

import sys
PYTHON_VERSION = sys.version_info.major

import gc
import os
import time
import logging
import itertools
import subprocess
from warnings import warn
from collections import defaultdict

try:
    if PYTHON_VERSION == 3:
        from edityper import toolkit
    elif PYTHON_VERSION == 2:
        import toolkit
        from itertools import imap as map
        range = xrange
    else:
        raise SystemExit("Please use Python 2.7 or 3.3 or higher")
except ImportError as error:
    raise SystemExit(error)


try:
    #   Use regex over re as it allows for fuzzy matching
    import regex
except ImportError as error:
    raise SystemExit(error)


SORT_ORDER = 'coordinate' # type: str
HD_HEADER = '@HD\tVN:1.5\tSO:' + SORT_ORDER # type: str
READ_GROUP_VALUES = { # type: Dict[str, str]
    #   Option      : header tag
    'read_center'   :   'CN',
    'read_library'  :   'LB',
    'read_platform' :   'PL',
    'read_sample'   :   'SM'
}

class SAM(object):
    """A SAM alignment
    This class contains nearly every piece of information that goes into a SAM alignment
    including (in __init__ order): query name (str), bitwise flag in integer form (int),
        reference name (str), position relative to reference (int), mapping quality (int)
        CIGAR string (str), next reference name (str), next position (int),
        template length (int), query sequence (str), query base quality (int)
    If providing position and pnext in 1-based coordinates, set 'zero_based' to 'False'

    This object will automatically find read group information

    Accessor methods are provided for some variables

    Three staticmethods are part of this class:
        validate_cigar:     make sure a CIGAR string is valid given a query sequence
        validate_flag:      make sure a integer flag is valid for SAM
        validate_quality:   ensure that the base quality is valid for a query sequence

    One class variable is provided: VALID_OPS is a set of valid CIGAR operations

    To print the SAM alignment like a SAM file, use str(sam)"""

    #   What are the valid CIGAR operation charactes?
    VALID_OPS = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'} # type: Set[str]

    #   What are valid tag types?
    _VALID_TAG_TYPES = {'A', 'B', 'H', 'Z', 'i', 'f'}

    @staticmethod
    def validate_cigar(cigar, sequence):
        """Validate a CIGAR string"""
        #   Broad check of CIGAR string
        valid_cigar = regex.compile(r'(^(\d+[%s])+$)' % ''.join(SAM.VALID_OPS)).search
        if not valid_cigar(cigar):
            raise ValueError("Invalid CIGAR string")
        #   Broad check CIGAR operations
        operations = ''.join((char for char in cigar if not char.isdigit())) # Get a string of the operations
        valid_ops = regex.compile(r'(^[%s]+$)' % ''.join(SAM.VALID_OPS)).search
        if not valid_ops(operations):
            raise ValueError("CIGAR string operations can only be from the following characters: '" + "', '".join(SAM.VALID_OPS + "'"))
        #   Validate hard clipping operations
        hard_indices = tuple(index for index, char in enumerate(operations) if char == 'H') # Get a list of indices where there was hard clipping
        if hard_indices: # This is only run if we had an 'H' in our operations
            try:
                if len(hard_indices) == 1:
                    assert hard_indices[0] == len(operations) - 1 or hard_indices[0] == 0
                elif len(hard_indices) == 2:
                    assert max(hard_indices) == len(operations) - 1 and min(hard_indices) == 0
                else: # Hard clipping can only happen up to two times. Something weird happened here
                    assert False
            except AssertionError:
                raise ValueError("Hard clipping can only be at the start or end of a sequence operation")
        #   Validate soft clipping operations
        if 'S' in operations:
            try:
                soft_clipping = regex.finditer(r'(S)', operations)
                soft_starts = tuple(soft.start() for soft in soft_clipping)
                #   If there is an S, any operation between it and the end of the CIGAR string must be an H
                start_ops = set(operations[min(soft_starts)::-1].replace('S', '')) # type: Set[str]
                final_ops = set(operations[max(soft_starts):].replace('S', '')) # type: Set[str]
                for end_ops in (start_ops, final_ops): # type: Set[str]
                    if end_ops:
                        assert len(end_ops) == 1 and 'H' in end_ops
            except AssertionError:
                raise ValueError("Between 'S' and the end of the CIGAR string, there can only be 'H' operations")
        #   Check to ensure that the number of bases in the CIGAR string is the same as the sequence length
        digits = regex.compile(r'(\d+)[MIS=X]').findall
        cigar_digits = map(int, digits(cigar)) # Get a list of the digits for the operations we need, as integers
        if sum(cigar_digits) != len(sequence):
            raise ValueError("The number of bases in the CIGAR string must equal the number of bases in the sequence")

    @staticmethod
    def validate_flag(flag):
        """Validate a SAM FLAG"""
        bits = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048} # type: Set[int]
        combinations = {seq for i in range(len(bits), 0, -1) for seq in itertools.combinations(bits, i) if sum(seq) == flag} # type: Set[int]
        if not combinations and flag != 0:
            raise ValueError("The FLAG must be the sum of any of the following: %s; not %s" % (', '.join(map(str, sorted(bits))), flag))

    @staticmethod
    def validate_quality(quality, sequence):
        """Validate the quality scores of the alignment"""
        try:
            if quality != '*':
                assert sequence != '*'
                assert len(quality) == len(sequence)
        except AssertionError:
            raise ValueError("If the quality is not a '*', the sequence must not be a '*' and the quality and sequence must be the same length")

    def __init__(
            self,
            qname='*', # type: str
            flag=4, # type: int
            rname='*', # type: str
            pos=0, # type: int
            mapq=255, # type: int
            cigar='*', # type: str
            rnext='*', # type: str
            pnext=0, # type: int
            tlen=0, # type: int
            seq='*', # type: str
            qual='*', # type: str
            zero_based=True # type: bool
    ):
        # type: (...) -> None
        try:
            #   Validate our inputs, namely the CIGAR string, FLAG, and QUALITY scores
            if cigar != '*':
                self.validate_cigar(cigar=cigar, sequence=seq)
            self.validate_flag(flag=flag)
            self.validate_quality(quality=qual, sequence=seq)
        except ValueError:
            raise
        self._qname = qname.split()[0] # toolkit.Read.name
        self._flag = flag # Bitwise flag in integer form
        self._rname = rname # From reference FASTA/Q file
        self._pos = pos + 1 if zero_based else pos # Were we given a zero-based coordinate? Add one for SAM specifications
        self._mapq = mapq # alignment.Alignment.score
        self._cigar = cigar # From write_sam.make_cigar_sequence()
        self._rnext = rnext # Not used for CRISPR program
        self._pnext = pnext + 1 if zero_based and rnext != '*' else pnext # Not used for CRISPR program
        self._tlen = tlen # Not used for CRISPR program
        self._seq = seq # From write_sam.make_cigar_sequence()
        self._qual = qual # toolkit.Read.qual
        #   Create a metadata set for holding optional tags, like RG and PG
        self._metadata = dict() # type: Dict[str, Tuple(str)]
        #   Make our read group
        self._set_read_group()

    def __repr__(self): # type: (None) -> str
        return self._qname + ':' + self._rname

    def __hash__(self): # type: (None) -> int
        return hash(self.__repr__())

    def __str__(self): # type: (None) -> str
        out = ( # type: Tuple[Any]
            self._qname,
            self._flag,
            self._rname,
            self._pos,
            self._mapq,
            self._cigar,
            self._rnext,
            self._pnext,
            self._tlen,
            self._seq,
            self._qual,
        )
        out = filter(lambda x: x is not None, out) # type: Iterable[Any]
        out = list(map(str, out)) # type: List[str]
        out += [':'.join((tag, tup[0], tup[1])) for tag, tup in sorted(self._metadata.items(), key=lambda tup: tup[0])]
        return '\t'.join(out)

    def __del__(self):
        del self._qname, self._flag, self._rname, self._pos, self._mapq, self._cigar, \
            self._rnext, self._pnext, self._tlen, self._seq, self._qual
        del self._metadata

    def __eq__(self, other):
        if not isinstance(other, SAM):
            return NotImplemented
        return self._rname == other._rname and self._pos == other._pos

    def __lt__(self, other):
        if not isinstance(other, SAM):
            return NotImplemented
        if self._rname == other._rname:
            return self._pos < other._pos
        else:
            return self._rname < other._rname

    def __le__(self, other):
        if not isinstance(other, SAM):
            return NotImplemented
        if self._rname == other._rname:
            return self._pos <= other._pos
        else:
            return self._rname <= other._rname

    def _set_read_group(self): # type: (None) -> None
        #   Specialized wrapper around SAM.add_metadata()
        if self._qname != '*':
            split_qname = self._qname.split(':')
            rgid = ':'.join((split_qname[0], split_qname[-4]))
            self.add_metadata(tag='RG', tag_type='Z', value=rgid)

    def _set_program(self, pg_id): # type: (str) -> None
        """Set this alignment's program ID
        Specialized wrapper around SAM.add_metadata()"""
        self.add_metadata(tag='PG', tag_type='Z', value=pg_id)

    def _get_lane(self): # type: (None) -> str
        """Get the lane for this alignment"""
        return ':'.join(self._qname.split(':')[:-3])

    def _get_read_group(self): # type: (None) -> str
        """Get the read group ID for this alignment
        Wrapper around SAM.get_metadata()"""
        read_group = self.get_metadata(tag='RG')
        if not read_group:
            warn("No read group set, returning 'None'")
        return read_group

    def _get_program(self): # type: (Optional[bool]) -> Union[str, None]
        """Get the program ID for this alignment
        Set 'full' to 'True' to get the program with 'PG:Z:' prefix
        Wrapper around SAM.get_metadata()"""
        return self.get_metadata(tag='PG')

    def add_metadata(self, tag, tag_type, value): # type: (str, str, Any) -> None
        """Set tag data"""
        #   Ensure we have a proper tag
        if tag_type not in self._VALID_TAG_TYPES:
            raise ValueError("'tag_type' must be one of '%s'" % "', '".join(self._VALID_TAG_TYPES))
        #   Type checking between tag_type and value
        if tag_type == 'A' and len(str(value)) != 1:
            raise ValueError("A tag of type 'A' must have a value that is a single character")
        elif tag_type == 'f' and not isinstance(value, (float, int)):
            raise ValueError("A tag of type 'f' must have a value of type 'float' or 'int'")
        elif tag_type == 'i' and not isinstance(value, int):
            raise ValueError("A tag of type 'i' must have a value of type 'int'")
        elif tag_type == 'Z' and not isinstance(value, str):
            raise ValueError("A tag of type 'Z' must have a value of type 'str'")
        #   Simple warning
        if tag in self._metadata:
            warn("A tag of '%s' was found in the metadata, overwriting" % tag)
        #   Add the metadata to our SAM object
        self._metadata[str(tag)] = (str(tag_type), str(value))

    def modify_read_group(self, modifier): # type: (str) -> Union[str, None]
        """Modify a read group ID to handle multiple lanes and avoid clashing"""
        try:
            #   Get the read group and modify it
            read_group = self.read_group # type: str
            new_rg_tag = read_group + str(modifier) # type: str
        except TypeError: # Concatenating a NoneType and str will raise a TypeError
            warn("Cannot modify an unset read group")
            new_read_group = None # type: NoneType
        else: # We properly modified our read group, now add it back to the metadata dictionary
            self.add_metadata(tag='RG', tag_type='Z', value=new_rg_tag)
            new_read_group = self.read_group # type: str
        return new_read_group

    def no_reference(self): # type: (None) -> None
        """Use this if there is no reference
        Will modify FLAG, RNAME, POS, MAPQ,
        CIGAR, RNEXT, and PNEXT fields
        to represent no reference found (basically unmapped)"""
        self._flag = 4 # type: int
        self._rname = '*' # type: str
        self._pos = 0 # type: int
        self._mapq = 255 # type: int
        self._cigar = '*' # type: str
        self._rnext = '*' # type: str
        self._pnext = 0 # type: int

    def get_rname(self): # type: (None) -> str
        """Get the reference feature name"""
        return self._rname

    def get_metadata(self, tag): # type: (str) -> Union[str, None]
        """Get metadata"""
        try:
            metadata = self._metadata.get(tag)[1]
        except TypeError:
            warn("Could not find a value in the metadata matching '%s', returning 'None'" % tag)
            metadata = None
        return metadata

    lane = property(fget=_get_lane, doc="Sample lane")
    read_group = property(fget=_get_read_group, doc="Sample read group")
    program = property(fget=_get_program, fset=_set_program, doc="Program used to make SAM file")


def get_genomic_location(bedfile): # type: (str) -> (str, int)
    """Get the chromosome and start location of the first line in a BED or GTF/GFF file, returns 1-based coordinates"""
    logging.info("Finding genomic location of the reference sequence")
    genomic_start = time.time() # type: float
    extension = os.path.splitext(bedfile)[1].lower() # type: str
    if extension in ('.gtf', '.gff'):
        logging.warning("Reading %s as GTF/GFF file", bedfile)
        start_index, position_modifier = 3, -1 # type: int, int
    elif extension == '.bed':
        logging.info("Reading %s", bedfile)
        start_index, position_modifier = 1, 0 # type: int, int
    else:
        logging.warning("Could not identify annotation file type, assuming BED format")
        extension = '.bed' # type: str
        start_index, position_modifier = 1, 0 # type: int, int
    with open(bedfile, 'r') as bfile: # type: _io.TextIOWrapper
        for line in bfile: # type: str
            line = line.strip() # type: str
            if line.startswith(('browser', 'track', '#')):
                continue
            line = line.split() # type: str
            try:
                chrom = line[0] # type: str
                genomic_location = int(line[start_index]) + position_modifier # type: int
            except (IndexError, ValueError):
                raise ValueError("Malformed %s file" % extension.upper().replace('.', ''))
            break
        else:
            raise ValueError("Could not parse %s file" % extension.upper().replace('.', ''))
    logging.debug("Finding the genomic location took %s seconds", round(time.time() - genomic_start, 3))
    return chrom, genomic_location


def make_sam_sequence(alignment, head=None, tail=None): # type: (str, Optional[int], Optional[int]) -> str
    """Make a SAM-ready aligned read"""
    if head and tail:
        sam_seq = alignment.read[head:tail + 1] # type: str
    else:
        sam_seq = alignment.read # type: str
    return sam_seq.replace('-', '')


def make_cigar(alignment, exact=False): # type: (al.Alignment, bool) -> (str, str)
    """Make a CIGAR string from an alignment and get
    a SAM-ready version of the aligned read"""
    #   Collect the reference and read sequences from the alignment
    reference = alignment.reference # type: str
    read = alignment.read # type: str
    if len(reference) != len(read):
        raise ValueError("The read and reference must be the same length")
    #   Get the start and end of our alignment, ignoring the '-'
    #   at the head and tail of the read as those
    #   are artifacts of alignment
    head, tail = toolkit.trim_interval(seq=read) # type: int, int
    cigar_ops = [] # type: List[str]
    #   Iterate through our range, comparing bases between read and reference
    for i in range(head, tail + 1): # type: int
        try:
            if reference[i] == '-' and read[i] == '-': # This shouldn't happen, but just in case
                continue
            #   Actual comparisons
            if reference[i] == '-': # Insertion relative to reference
                cigar_ops.append('I')
            elif read[i] == '-': # Deletion relative to reference
                cigar_ops.append('D')
            elif exact and read[i] == reference[i]: # Sequence match, when exact CIGAR is needed
                cigar_ops.append('=')
            elif exact and read[i] != reference[i]: # Sequence mismatch, when exact CIGAR is needed
                cigar_ops.append('X')
            else: # Alignment match (not base match)
                cigar_ops.append('M')
        except IndexError:
            break
    #   Now, turn our list of operations into a proper cigar string
    counter = 0 # type: int
    cigar_string = '' # type: str
    prev_op = cigar_ops[0]
    for index, operation in enumerate(cigar_ops):
        #   Are we still on the same operation?
        if operation == prev_op:
            counter += 1 # If so, increment our counter
            if index != len(cigar_ops) - 1: # If this isn't the last operation, iterate
                continue
        #   Go straight here if new operation
        #   or come here if this is the last operation
        #   and the previous operation continued
        cigar_string += str(counter) + prev_op
        counter = 1 # type: int
        #   If this is the last operation AND not a continuation of the previous operation
        if index == len(cigar_ops) - 1 and operation != prev_op:
            #   Have to add it to the end of our cigar string
            cigar_string += str(counter) + operation
        prev_op = operation # type: str
    #   Validate our cigar string and get a SAM-ready version of our read
    SAM.validate_cigar(cigar=cigar_string, sequence=make_sam_sequence(alignment=alignment, head=head, tail=tail))
    return cigar_string


def make_read_group(sam_lines, conf_dict): # type: (List[SAM], Dict[Any]) -> Tuple[str]
    """Make read group @RG headers"""
    #   Two holding dictionaries,
    #   one for all IDs and PUs and
    #   one for unique ID/PU pairs
    raw_ids_pus = defaultdict(set) # type: defaultdict[Set[str]]
    rg_ids_pus = dict() # type: Dict[str, str]
    for sam in sam_lines:
        raw_ids_pus[sam.read_group].add(sam.lane)
    rg_lines = list() # type: List[str]
    for rgid, rgpus in raw_ids_pus.items():
        #   If there are more than one PU for this ID
        if len(rgpus) > 1:
            logging.error("There were multiple lanes found for read group ID: %s", rgid)
            #   Modify with a simple _#
            mod_counter = 1 # type: int
            #   For every PU for this ID
            for rgpu in rgpus: # type: str
                #   Find SAM lines with the same PU
                for sam in sam_lines: # type: SAM
                    if sam.lane == rgpu:
                        #   Modify the RGID for this SAM
                        new_rgid = sam.modify_read_group(modifier='_' + str(mod_counter)) # type: str
                #   Increment our mod counter for the next PU
                #   and add the ID/PU pair to our final dictionary
                mod_counter += 1
                rg_ids_pus[new_rgid] = rgpu
        else:
            #   Otherwise, just roll with what was found
            rg_ids_pus[rgid] = rgpus.pop()
    #   For every final (including modified) ID/PU pair we have, make an @RG header for it
    for rgid, rgpu in rg_ids_pus.items():
        line = '@RG\tID:' + rgid + '\tPU:' + rgpu # type: str
        rg_lines.append(line)
    header_tags = [] # type: List[str]
    for option, tag in READ_GROUP_VALUES.items(): # type: str, str
        if option in conf_dict:
            complete_tag = ':'.join((tag, conf_dict[option])) # type: str
            header_tags.append(complete_tag)
    if header_tags:
        for line in rg_lines: # type: str
            line += '\t' + '\t'.join(header_tags)
    return tuple(rg_lines)


def make_sequence_header(
        sam_lines, # type: Iterable[SAM]
        ref_seq_dict # type: Dict[str, Tuple[str, int]]
):
    # type: (...) -> Tuple[str]
    """Make reference @SQ headers
    sam_lines       :   Iterable[SAM objects]
    ref_seq_dict    :   Dict[ref seq name, Tuple[ref seq, global modifier (usually 0)]]"""
    rnames = {sam.get_rname() for sam in sam_lines} # type: Set[str]
    for rname in rnames: # type: str
        #   If the RNAME is not in our reference,
        #   modify any SAM lines that map to it as unmapped
        #   and remove the RNAME for our set of RNAMEs so
        #   that we do not make a @SQ header for it
        if rname == '*':
            continue
        if rname not in ref_seq_dict:
            logging.error("Could not find reference sequence %s", rname)
            for sam in sam_lines: # type: SAM
                if sam.get_rname() == rname:
                    sam.no_reference()
            rnames.discard(rname)
    #   Make @SQ headers for the RNAMEs we found in our reference
    sq_header = [] # type: List[str]
    for rname in sorted(rnames): # type: str
        if rname == '*':
            continue
        seq, mod = ref_seq_dict[rname] # type: str, int
        # seq = ref_seq_dict[rname].replace('-', '') # type: str
        seq = seq.replace('-', '') # type: str
        line = ('@SQ', 'SN:' + rname, 'LN:%d' % (len(seq) + mod)) # type: generator
        sq_header.append('\t'.join(line))
    return tuple(sq_header)


def calc_read_pos(alignment, genomic_start=0): # type: (alignment.Alignment, int) -> int
    """Calculate read position"""
    head, _ = toolkit.trim_interval(seq=alignment.read) # type: int, int
    cigar = make_cigar(alignment=alignment) # type: str
    first_consumed = regex.search(r'(\d+[MDN=X])', cigar) # type: _regex.Match
    if not first_consumed:
        raise ValueError("Read did not align to reference")
    cigar = cigar[:first_consumed.start()] # type: str
    not_consumed = regex.findall(r'\d+', cigar) # type: List[str]
    not_consumed = map(int, not_consumed) # type: Iterable[str]
    # import code; code.interact(local=locals()); sys.exit()
    return head + genomic_start + sum(not_consumed)


def bam(fastq_name, samfile, samtools, index_type): # type: (str, str, str, str) -> None
    """Use SAMtools to convert a SAM to BAM file"""
    logging.info("FASTQ %s: Converting SAM to BAM", fastq_name)
    bam_start = time.time() # type: float
    #   Make BAM name
    bamfile = samfile.replace('sam', 'bam') # type: str
    view_cmd = (samtools, 'view -bhS', samfile, '>', bamfile) # type: Tuple[str]
    index_arg = 'c' if index_type == 'csi' else 'b' # type: str
    index_cmd = '%(samtools)s index -%(arg)s %(bamfile)s' % {'samtools': samtools, 'arg': index_arg, 'bamfile': bamfile} # type: str
    logging.info("FASTQ %s: Writing BAM to %s", fastq_name, bamfile)
    subprocess.call(' '.join(view_cmd), shell=True)
    gc.collect()
    logging.info("FASTQ %s: Indexing BAM file", fastq_name)
    logging.debug("FASTQ %s: Making %s indices", fastq_name, index_type)
    subprocess.call(index_cmd, shell=True)
    gc.collect()
    logging.debug("FASTQ %s: Converting SAM to BAM took %s seconds", fastq_name, round(time.time() - bam_start, 3))
    logging.debug("FASTQ %s: Removing SAM file, leaving only BAM file", fastq_name)
    os.remove(samfile)
