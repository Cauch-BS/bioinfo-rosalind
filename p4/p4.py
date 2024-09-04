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
    seq = np.frombuffer(seq_file.readline().strip().encode().translate(nuc_to_inc), dtype=np.uint8)
    windowlen = int(seq_file.readline().strip())
    seqlen = seq.shape[0]
    newshape = (seqlen - windowlen + 1, windowlen)
    newstrides = (seq.strides[0], seq.strides[0])
    windowed_seqs = as_strided(seq, shape=newshape, strides=newstrides)
    unique_rows, counts = np.unique(windowed_seqs, axis=0, return_counts=True)
    max_count_nuc_idx = nucleotides[unique_rows[counts == np.max(counts)]]
    max_count_nucs = np.array([''.join(row) for row in max_count_nuc_idx])
    print(*max_count_nucs)