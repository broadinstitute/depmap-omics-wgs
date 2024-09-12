from math import floor

import pysam


def calc_gc_percentage(
    chrom: str,
    pos: int,
    ref: str,
    variant_class: str,
    window_size: int,
    fasta_handle: pysam.FastaFile,
) -> float:
    if variant_class == "deletion":
        # shift to the right to account for ref/alt alleles containing left bp context
        pos += 1
    elif variant_class == "substitution" and len(ref) > 2:
        # center window around middle of sequence (e.g. move 1 to the right for
        # sequences of length 3 or 4)
        pos += floor((len(ref) - 1) / 2)

    # calculate the start and end of a symmetric window around pos
    half_window = window_size // 2
    start = max(1, pos - half_window)
    end = pos + half_window

    # extract the sequence from the fasta file
    window_seq = fasta_handle.fetch(chrom, start - 1, end)

    # calculate GC percentage
    return (window_seq.count("G") + window_seq.count("C")) / len(window_seq)
