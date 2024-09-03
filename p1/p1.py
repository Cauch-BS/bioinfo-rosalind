from dotenv import dotenv_values
import numpy as np

# Load the DNA sequence from the file
datadir = dotenv_values(".env")["DATA_DIR"]
loadfile = "/rosalind_dna.txt"
nucseq = open(f"{datadir}{loadfile}", "r", encoding=None)
nuc_to_inc = bytes.maketrans(b"ACGT", b"\x00\x01\x02\x03")

with nucseq as seq_file:
    seq = np.eye(4, dtype = np.int32)[
        np.frombuffer(seq_file.read().strip().encode().translate(nuc_to_inc), dtype=np.uint8)
    ].sum(axis=0)
    print(*seq)
    
    
