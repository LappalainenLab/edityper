#!/usr/bin/env python2

'''Alignment functions for the CRISPR program'''

from __future__ import print_function

import sys
if sys.version_info.major is not 2 and sys.version_info.minor is not 7:
    sys.exit("Please use Python 2.7 for this module: " + __name__)


from collections import defaultdict, namedtuple

try:
    import NW_py as NW
    import genetic_toolpack as toolpack
except ImportError as error:
    sys.exit("Please keep this program in it's directory to load custom modules: " + error.message)


ReadSummary = namedtuple("ReadSummary", ('sequence', 'reads'))

class Alignment(object):

    """An alignment just to make my life easier"""

    @classmethod
    def from_alignment_dict(cls, align_dict): # type: (Dict[str, Any]) -> Alignment
        """Create an Alignment from another Alignment's __dict__ method"""
        #   Use __init__ for parts of the new alignment
        #   As these are required arguments, just call them from the dictionary
        self = cls( # type: Alignment
            ref_align=align_dict['_ref'],
            read_align=align_dict['_read'],
            score=align_dict['_score'],
            names=align_dict['_names']
        )
        #   Get the rest of the information
        #   These are not required, so use dict.get() to
        #   ask if the key is in our alignment dictionary
        #   Set the value in our instance if so, otherwise
        #   set to None (None is default value of dict.get())
        self._sources = align_dict.get('_sources') # type: Optional[str]
        self._unaligned = align_dict.get('_unaligned') # type: Optional[str]
        self._num_reads = align_dict.get('_num_reads') # type: Optional[int]
        self._nmdel = align_dict.get('_nmdel') # type: Optional[int]
        self._nmins = align_dict.get('_nmins') # type: Optional[int]
        self._nmmis = align_dict.get('_nmmis') # type: Optional[int]
        return self

    @staticmethod
    def from_dict(other_dict): # type: (Dict[str, Any]) -> Alignment
        """Create an Alignment from a dictionary with the following keys:
        Required: 'ref', 'read', 'score', 'names'
        Optional: 'unaligned', 'num_reads', 'nmdel', 'nmins', 'nmmis'"""
        try:
            return Alignment.from_alignment_dict(align_dict={'_' + key: value for key, value in other_dict.items()})
        except TypeError:
            raise KeyError("An Alignment only has the following attributes: 'ref', 'read', 'score', 'names', 'sources' 'unaligned', 'num_reads', 'nmdel', nmins', and 'nmmis'")

    def __init__(
            self,
            ref_align, # type: str
            read_align, # type: str
            score, # type: str
            names=None # type: Optional[Iterable[str]]
    ):
        # type: (...) -> None
        self._ref = ref_align # type: str
        self._read = read_align # type: str
        self._score = score # type: str
        self._names = names # type: Iterable[str]
        self._sources = None # type: Optional[Iterable[str]]
        self._unaligned = None # type: Optional[str]
        self._num_reads = None # type: Optional[int]
        self._nmdel = None # type: Optional[int]
        self._nmins = None # type: Optional[int]
        self._nmmis = None # type: Optional[int]

    def __repr__(self): # type: (None) -> str
        return self._unaligned

    def set_unaligned(self, sequence): # type: (str) -> None
        """Provide an unaligned sequence"""
        self._unaligned = sequence

    def set_sources(self, sources): # type: Iterable(str) -> None
        """Set the source FASTQ files for this alignment"""
        self._sources = {source for source in sources}

    def set_stats(self, num_reads, nm_del, nm_ins, nm_mis):  # type: (int, int, int, int) -> None
        """Set the stats"""
        self._num_reads = num_reads
        self._nmdel = nm_del
        self._nmins = nm_ins
        self._nmmis = nm_mis

    def get_aligned_reference(self): # type: (None) -> str
        """Get the aligned reference"""
        return self._ref

    def get_aligned_read(self): # type: (None) -> str
        """Get the aligned read"""
        return self._read

    def get_score(self): # type: (None) -> int
        """Get the alignment score"""
        return self._score

    def get_names(self): # type: (None) -> Tuple[str]
        """Get the name of this alignment"""
        return tuple(self._names)

    def get_sources(self): # type: (None) -> Tuple[str]
        """Get the FASTQ files that support this alignment"""
        return tuple(self._sources)

    def get_stats(self): # type: (None) -> (int, int, int, int)
        """Get the stats for this alignment"""
        return self._num_reads, self._nmdel, self._nmins, self._nmmis

    def get_unaligned(self): # type: (None) -> str
        """Get the unaligned sequence"""
        return self._unaligned


def sort_reads_by_length(reads_dict): # type: (Dict[str, List[toolpack.Read]]) -> Dict[int, List[ReadSummary]]
    """Sort a list of reads by their length"""
    reads_by_length = defaultdict(list)
    for seq, reads_list in reads_dict.items(): # type: str, List[toolpack.Read]
        length = len(seq) # type: int
        seq_info = ReadSummary(sequence=seq, reads=tuple(read for read in reads_list)) # type: ReadSummary
        reads_by_length[length].append(seq_info)
        reads_by_length[length].sort(key=lambda summary: summary.sequence)
    return reads_by_length


def align_recurse(
        reads_by_length, # type: Dict[int, List[ReadSummary]]
        reference, # type: str
        gap_open, # type: int
        gap_ext, # type: int
):
    # type: (...) -> Dict[int, List[Alignment]]
    """Align using the recurisve method"""
    #   Keep track of matrix re-used lines
    alignments = defaultdict(list) # type: defaultdict[List]
    reuse = 0 # type: int
    for length, reads_list in reads_by_length.items(): # type: int, List[ReadSummary]
        count, temp = 0, '' # type: int, str
        total = len(reads_list)
        for summary in reads_list: # type: ReadSummary
            count += 1
            seq = summary.sequence # type: str
            names = tuple(read.get_readid() for read in summary.reads) # type: Tuple[str]
            fastqs = {read.get_source() for read in summary.reads} # type: Set[str]
            if not temp: # basically, first sequence to be aligned
                al_ref, al_read, score = NW.align_aff_mem(reference, seq, gap_open, gap_ext, 0, 0) # type: str, str, int
                aligned = Alignment(ref_align=al_ref, read_align=al_read, score=score, names=names) # type: Alignment
                aligned.set_unaligned(sequence=seq)
                aligned.set_sources(sources=fastqs)
                alignments[length].append(aligned)
                temp = seq # type: str
                continue
            index = toolpack.sim_seq(seq1=al_ref, seq2=al_read) # type: int
            reuse += index
            if count is total:
                al_ref, al_read, score = NW.align_aff_mem(reference, seq, gap_open, gap_ext, index, 1) # type: str, str, int
            else:
                al_ref, al_read, score = NW.align_aff_mem(reference, seq, gap_open, gap_ext, index, 0) # type: str, str, int
            aligned = Alignment(ref_align=al_ref, read_align=al_read, score=score, names=names)
            aligned.set_unaligned(sequence=seq)
            aligned.set_sources(sources=fastqs)
            alignments[length].append(aligned)
            temp = seq # type: str
    return dict(alignments)


def alignment_by_fastq(
        alignments, # type: Iterable[alignment.Alignment]
):
    # type: (...) -> Dict[str, List[alignment.Alignment]]
    """Organize alignments by the FASTQ files that contained the reads supporting each alignment"""
    sorted_alignments = defaultdict(list)
    for aligned in alignments:
        for source in aligned.get_sources():
            sorted_alignments[source].append(aligned)
    return sorted_alignments
