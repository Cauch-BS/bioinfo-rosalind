from dotenv import dotenv_values

# Load the DNA sequence from the file
datadir = dotenv_values(".env")["DATA_DIR"]
loadfile = "/rosalind_rna.txt"
with open(f"{datadir}{loadfile}", "r", encoding=None) as seqfile:
    seq = seqfile.read().strip().replace("T", "U")
    print(seq)
