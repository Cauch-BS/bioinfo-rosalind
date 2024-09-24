from dotenv import dotenv_values
import numpy as np
from numpy.lib.stride_tricks import as_strided

# Load the DNA sequence from the file
datadir = dotenv_values(".env")["DATA_DIR"]
loadfile = "/rosalind_ba1b.txt"
nucseq = open(f"{datadir}{loadfile}", "r", encoding=None)
nuc_to_inc = bytes.maketrans(b"ACGT", b"\x00\x01\x02\x03")
nucleotides = np.array(list("ACGT"), dtype = 'U1')

with nucseq as seq_file:
    motif_length, max_mismatch = map(int, seq_file.readline().strip().split())
    seq_list: list[np.ndarray] = []
    for line in seq_file:
        seq = np.frombuffer(line.strip().encode().translate(nuc_to_inc), dtype=np.uint8)
        windowlen = motif_length
        seqlen = seq.shape[0]
        newshape = (seqlen - windowlen + 1, windowlen)
        newstrides = (seq.strides[0], seq.strides[0])
        windowed_seqs = as_strided(seq, shape=newshape, strides=newstrides)
        seq_list.append(windowed_seqs)
    unique_motifs, _ = np.unique(np.concatenate(seq_list), axis=0, return_counts=False)
    motif_count = np.full(unique_motifs.shape[0], False, dtype=bool)
    for seq in seq_list:
        motif_count = np.sum(unique_motifs != seq, axis = 1)
        if np.min(motif_count) <= max_mismatch:
            motif_count[np.argmin(motif_count)] = True
