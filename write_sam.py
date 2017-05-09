#!/usr/bin/env python2

"""Create a SAM file from alignment data for the CRISPR program"""

from __future__ import print_function, division

import sys
if sys.version_info.major is not 2 and sys.version_info.minor is not 7:
    sys.exit("Please use Python 2.7 for this module: " + __name__)


import logging
import itertools
from copy import copy
from warnings import warn
from collections import defaultdict

try:
    import genetic_toolpack as toolpack
except ImportError as error:
    sys.exit("Please keep this program in it's directory to load custom modules: " + error.message)


try:
    #   Use regex over re as it allows for fuzzy matching
    import regex
except ImportError as error:
    sys.exit("Please install 'regex' for this module: " + error.message)


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
        hard_indecies = tuple(index for index, char in enumerate(operations) if char == 'H') # Get a list of indecies where there was hard clipping
        if hard_indecies: # This is only run if we had an 'H' in our operations
            try:
                if len(hard_indecies) == 1:
                    assert hard_indecies[0] == len(operations) - 1 or hard_indecies[0] == 0
                elif len(hard_indecies) == 2:
                    assert max(hard_indecies) == len(operations) - 1 and min(hard_indecies) == 0
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
        combinations = {seq for i in xrange(len(bits), 0, -1) for seq in itertools.combinations(bits, i) if sum(seq) == flag} # type: Set[int]
        if not combinations:
            raise ValueError("The FLAG must be the sum of any of the following: " + ', '.join(map(str, sorted(bits))))

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
        self._qname = qname # genetic_toolpack.Read._readid
        self._flag = flag # Bitwise flag in integer form
        self._rname = rname # From reference FASTA/Q file
        self._pos = pos + 1 if zero_based else pos # Were we given a zero-based coordinate? Add one for SAM specifications
        self._mapq = mapq # alignment.Alignment._score
        self._cigar = cigar # From write_sam.make_cigar_sequence()
        self._rnext = rnext # Not used for CRISPR program
        self._pnext = pnext + 1 if zero_based else pnext # Not used for CRISPR program
        self._tlen = tlen # Not used for CRISPR program
        self._seq = seq # From write_sam.make_cigar_sequence()
        self._qual = qual # genetic_toolpack.Read._quality
        #   Create a metadata set for holding optional tags, like RG and PG
        # self._metadata = set() # type: Set[str]
        self._metadata = dict() # type: Dict[str, Tuple(str)]
        #   Make our read group
        self._set_read_group()

    def __repr__(self): # type: (None) -> str
        return self._qname + ':' + self._rname

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
        out = filter(lambda x: x is not None, out) # type: Tuple[Any]
        out = map(str, out) # type: List[str]
        out += [':'.join((tag, tup[0], tup[1])) for tag, tup in sorted(self._metadata.items(), key=lambda tup: tup[0])]
        return '\t'.join(out)

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

    def set_program(self, pg_id): # type: (str) -> None
        """Set this alignment's program ID
        Specialized wrapper around SAM.add_metadata()"""
        self.add_metadata(tag='PG', tag_type='Z', value=pg_id)

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
            read_group = self.get_read_group() # type: str
            new_rg_tag = read_group + str(modifier) # type: str
        except TypeError: # Concatenating a NoneType and str will raise a TypeError
            warn("Cannot modify an unset read group")
            new_read_group = None # type: NoneType
        else: # We properly modified our read group, now add it back to the metadata dictionary
            self.add_metadata(tag='RG', tag_type='Z', value=new_rg_tag)
            new_read_group = self.get_read_group() # type: str
        finally:
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
            metadata = tuple(meta for meta in self._metadata if tag in meta)[0] # type: str
        except IndexError:
            warn("Could not find a value in the metadata matching '%s', returning 'None'" % tag)
            metadata = None
        finally:
            return metadata

    def get_lane(self): # type: (None) -> str
        """Get the lane for this alignment"""
        return ':'.join(self._qname.split(':')[:-3])

    def get_read_group(self): # type: (None) -> str
        """Get the read group ID for this alignment
        Wrapper around SAM.get_metadata()"""
        read_group = self.get_metadata(tag='RG')
        try:
            read_group = ':'.join(read_group.split(':')[2:])
        except AttributeError: # NoneTypes don't have a split attribute
            warn("No read group set, returning 'None'")
        return read_group

    def get_program(self, full=False): # type: (Optional[bool]) -> Union[str, None]
        """Get the program ID for this alignment
        Set 'full' to 'True' to get the program with 'PG:Z:' prefix
        Wrapper around SAM.get_metadata()"""
        program = self.get_metadata(tag='PG')
        if not full:
            try: program = ':'.join(program.split(':')[2:])
            except AttributeError: pass
        return program


def make_sam_sequence(alignment, head=None, tail=None): # type: (str, Optional[int], Optional[int]) -> str
    """Make a SAM-ready aligned read"""
    if head and tail:
        sam_seq = alignment.get_aligned_read()[head:tail + 1] # type: str
    else:
        sam_seq = alignment.get_aligned_read() # type: str
    return sam_seq.replace('-', '')


def make_cigar(alignment, exact=False): # type: (al.Alignment, bool) -> (str, str)
    """Make a CIGAR string from an alignment and get
    a SAM-ready version of the aligned read"""
    #   Collect the reference and read sequences from the alignment
    reference = alignment.get_aligned_reference() # type: str
    read = alignment.get_aligned_read() # type: str
    if len(reference) != len(read):
        raise ValueError("The read and reference must be the same length")
    #   Get the start and end of our alignment, ignoring the '-'
    #   at the head and tail of the read as those
    #   are artifacts of alignment
    head, tail = toolpack.trim_interval(seq=read) # type: int, int
    cigar_ops = [] # type: List[str]
    #   Iterate through our range, comparing bases between read and reference
    for i in xrange(head, tail + 1): # type: int
        if reference[i] == '-' and read[i] == '-': continue # This shouldn't happen, but just in case
        #   Actual comparisons
        if reference[i] == '-': cigar_ops.append('I') # Insertion relative to reference
        elif read[i] == '-': cigar_ops.append('D') # Deletion relative to reference
        elif exact and read[i] == reference[i]: cigar_ops.append('=') # Sequence match, when exact CIGAR is needed
        elif exact and read[i] != reference[i]: cigar_ops.append('X') # Sequence mismatch, when exact CIGAR is needed
        else: cigar_ops.append('M') # Alignment match (not base match)
    #   Now, turn our list of operations into a proper cigar string
    counter = 0 # type: int
    cigar_string = '' # type: str
    prev_op = cigar_ops[0]
    for index, operation in enumerate(cigar_ops):
        #   Are we still on the same operation?
        if operation == prev_op:
            counter += 1 # If so, increment our counter
            if index != len(cigar_ops) - 1: continue # If this isn't the last operation, iterate
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


def calc_read_pos(alignment): # type: alignment.Alignment -> int
    """Calculate alignment read position"""
    cigar = make_cigar(alignment=alignment, exact=True) # type: str
    reference = alignment.get_aligned_reference() # type: str
    read = alignment.get_aligned_read() # type: str
    if len(reference) != len(read):
        raise ValueError("The read and reference must be the same length")
    reference = reference.replace('-', '') # type: str
    read = read.replace('-', '') # type: str
    SAM.validate_cigar(cigar=cigar, sequence=read)
    #   Split our cigar string into a list of component operations
    operations = regex.findall(r'(\d+[%s])' % ''.join(SAM.VALID_OPS), cigar) # type Tuple[str]
    #   Create a regex for finding digits
    digits = regex.compile(r'(\d+)[%s]' % ''.join(SAM.VALID_OPS)).findall
    #   Start our pattern for the final regex search
    pattern = '(' # type: str
    #   Create this iterable to validate the checks below
    possible_references = tuple(itertools.repeat(reference, 2)) # type: Tuple[str]
    #   Count the number of types we modify the reference
    #   Use this for fuzzy matching bounds
    #   when searching for a new reference
    num_ref_adjust = 0 # type: int
    #   Start parsing the CIGAR string operations
    for oper in operations: # type: str
        this_op = ''.join((char for char in oper if not char.isdigit())) # type: str
        this_digit = int(digits(oper)[0]) # type: int
        #   Pattern/read/reference adjustment
        #   If sequence match
        if this_op in {'='}:
            #   Add to the pattern the sequence that matches
            #   Adjust the read because we've covered this sequence
            pattern += read[:this_digit]
            read = read[this_digit:] # type: str
        #   If sequence mismatch, but alignment match
        elif this_op in {'X'}:
            #   Start a counter so that we only go up to the number of sequence mismatches
            i = 0 # type: int
            while i < this_digit:
                #   Add to the pattern the base that doesn't match
                #   Remove this base from the read
                #   Increment the counter
                pattern += '[^%s]' % read[0]
                read = read[1:] # type: str
                i += 1
        #   If we have part of our read that isn't in our reference
        elif this_op in {'I', 'S', 'P'}: # Insertion in read relative to reference, soft clipping of read, or silent deletion of reference
            #   Remove these bases from the read and continue, don't check
            read = read[this_digit:] # type: str
            continue
        elif this_op in {'D'}: # Deletion in read relative to reference
            adjustments = list(regex.finditer(r'%s.{%s})' % (pattern, this_digit), ref) for ref in possible_references) # type: List
            adjustments = (match.groups() for adj in adjustments for match in adj) # type: generator
            adjustments = itertools.chain.from_iterable(adjustments) # type: itertools.chain
            adj_ref = copy(possible_references) # type: List[str]
            for adj in adjustments: # type: str
                for ref in adj_ref: # type: str
                    temp = list() # type: List[str]
                    try:
                        ref.index(adj)
                        ref.replace(adj, adj[:-this_digit])
                        temp.append(ref)
                    except ValueError:
                        continue
                adj_ref = copy(temp) # type: List[str]
            possible_references = copy(adj_ref) # type: List[str]
            num_ref_adjust += this_digit
        else: # Hard clipping, sequence not found in seq
            pass
        possible_references = (regex.findall(r'%s.*)' % pattern, ref) for ref in possible_references) # type: generator
        possible_references = tuple(itertools.chain.from_iterable(possible_references)) # type: Tuple[str]
        if len(possible_references) == 1:
            break
        elif len(possible_references) == 0:
            raise ValueError("Something happened")
    if len(possible_references) > 1:
        logging.error("Multiple mapping is currently not supported, setting position to '-8' for fixing later")
        return -8
    adjusted_ref = possible_references[0] # type: str
    return regex.search(r'(%s){e<=%d}' % (adjusted_ref, num_ref_adjust), reference, regex.BESTMATCH).start()


def make_read_group(sam_lines, conf_dict): # type: (List[SAM], Dict[Any]) -> Tuple[str]
    """Make read group @RG headers"""
    #   Two holding dictionaries,
    #   one for all IDs and PUs and
    #   one for unique ID/PU pairs
    raw_ids_pus = defaultdict(set) # type: defaultdict[Set[str]]
    rg_ids_pus = dict() # type: Dict[str, str]
    for sam in sam_lines:
        raw_ids_pus[sam.get_read_group()].add(sam.get_lane())
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
                    if sam.get_lane() == rgpu:
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
        ref_seq_dict # type: Dict[str, str]
):
    """Make reference @SQ headers"""
    rnames = {sam.get_rname() for sam in sam_lines} # type: Set[str]
    for rname in rnames: # type: str
        #   If the RNAME is not in our reference,
        #   modify any SAM lines that map to it as unmapped
        #   and remove the RNAME for our set of RNAMEs so
        #   that we do not make a @SQ header for it
        if rname not in ref_seq_dict:
            logging.error("Could not find reference sequence %s", rname)
            for sam in sam_lines: # type: SAM
                if sam.get_rname() == rname:
                    sam.no_reference()
            rnames.discard(rname)
    #   Make @SQ headers for the RNAMEs we found in our reference
    sq_header = [] # type: List[str]
    for rname in sorted(rnames): # type: str
        seq = ref_seq_dict[rname].replace('-', '') # type: str
        line = ('@SQ', 'SN:' + rname, 'LN:%d' % len(seq)) # type: generator
        sq_header.append('\t'.join(line))
    return tuple(sq_header)


def make_sam(
        alignment, # type: alignment.Alignment
        sam_seq, # type: str
        cigar, # type: str
        bit_flag, # type: int
        ref_name, # type: str
        position, # type: int
        reads_dict # type: Dict[str, List[toolpack.Read]]
):
    # type: (...) -> Tuple[SAM]
    """Make a collection of SAM objects"""
    sam_lines = list() # type: List[str]
    try:
        reads = {read.get_readid(): read for read in reads_dict[alignment.get_unaligned()]} # type: Dict[str, toolpack.Read]
    except KeyError:
        logging.error(
            "Cannot find the reads that support the alignment for %s, some data will be missing",
            alignment.get_unaligned()
        )
        reads = dict() # Used so that only one other try/except clause is needed, not two.
    read_names = sorted((name for name in alignment.get_names())) # type: List[str]
    for readid in read_names: # type: str
        try:
            quality = reads[readid].get_quality() # type: str
        except KeyError:
            quality = '*' # type: str
        sam = SAM(
            qname=readid.replace('@', '').split()[0],
            flag=bit_flag,
            rname=ref_name,
            pos=position,
            # mapq=alignment.get_score(),
            mapq=255,
            cigar=cigar,
            rnext='*',
            pnext=0,
            tlen=0,
            seq=sam_seq,
            qual=quality
        )
        sam.add_metadata(tag='AS', tag_type='i', value=alignment.get_score())
        sam_lines.append(sam)
    #   Add unaligned reads
    #   This can be here because unaligned will return empty
    #   if reads is just a generic dictionary
    unaligned = {readid: read for readid, read in reads.items() if readid not in read_names} # type: Dict[str, toolpack.Read]
    if unaligned:
        logging.error("Found %s unaligned reads", len(unaligned))
        for readid, read in unaligned.items(): # type: str, toolpack.Read
            sam = SAM(
                qname=readid.replace('@', '').split()[0],
                flag=4,
                seq=read.get_sequence(),
                qual=read.get_quality()
            )
            sam_lines.append(sam)
    return tuple(sam_lines)
