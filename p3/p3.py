from dotenv import dotenv_values
import numpy as np
from numpy.lib.stride_tricks import as_strided

# Load the DNA sequence from the file
datadir = dotenv_values(".env")["DATA_DIR"]
loadfile = "/rosalind_ba1a.txt"
nucseq = open(f"{datadir}{loadfile}", "r", encoding=None)
nuc_to_inc = bytes.maketrans(b"ACGT", b"\x00\x01\x02\x03")

with nucseq as seq_file:
    seq = np.frombuffer(seq_file.readline().strip().encode().translate(nuc_to_inc), dtype=np.uint8)
    window = np.frombuffer(seq_file.readline().strip().encode().translate(nuc_to_inc), dtype = np.uint8)
    seqlen = seq.shape[0]
    windowlen = window.shape[0]
    newshape = (seqlen - windowlen + 1, windowlen)
    newstrides = (seq.strides[0], seq.strides[0])
    windowed_seqs = as_strided(seq, shape=newshape, strides=newstrides)
    print(np.sum(np.all(windowed_seqs == window, axis=1)))