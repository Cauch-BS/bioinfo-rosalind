from itertools import product

from dotenv import dotenv_values
import numpy as np
from numpy.lib.stride_tricks import as_strided

def check_motif(strided_seq: np.ndarray, 
                motif_candidates: np.ndarray,
                max_mismatch: int) -> np.ndarray:
    mismatches = np.sum(strided_seq[:, np.newaxis] != motif_candidates, axis = 2)
    mask = np.any(mismatches <= max_mismatch, axis = 0)
    pruned_motif_candidates = motif_candidates[mask]

    # return the pruned motif candidates
    return pruned_motif_candidates

def main(mode: str):
    # Load the DNA sequence from the file
    datadir = dotenv_values(".env")["DATA_DIR"]
    loadfile = "rosalind_ba2a.txt"
    nucseq = open(f"{datadir}{loadfile}", "r", encoding=None)
    nuc_to_inc = bytes.maketrans(b"ACGT", b"\x00\x01\x02\x03")
    nucleotides = np.array(list("ACGT"), dtype = 'U1')

    with nucseq as seq_file:
        motif_length, max_mismatch = map(int, seq_file.readline().strip().split())
        strided_seqs = list()
        if 'debug' in mode:
            seq_list = seq_file.readline().strip().split()
        elif 'actual' in mode:
            seq_list = []
            for line in seq_file:
                seq_list.append(line.strip())
        else:
            raise ValueError(f"Invalid mode: {mode}")
        
        
        for line in seq_list:
            seq = np.frombuffer(line.strip().encode().translate(nuc_to_inc), dtype=np.uint8)
            seqlen = seq.shape[0]
            newshape = (seqlen - motif_length + 1, motif_length)
            newstrides = (seq.strides[0], seq.strides[0])
            strided_seq = as_strided(seq, shape=newshape, strides=newstrides)
            strided_seqs.append(strided_seq)

        motif_candidates = np.array(list(product(range(4), repeat=motif_length)))

        for strided_seq in strided_seqs:
            motif_candidates = check_motif(
                strided_seq,
                motif_candidates,
                max_mismatch
            )
        motifs = nucleotides[motif_candidates]
        motif_list = [
            ''.join(row) for row in motifs
        ]
        print(*motif_list)

if __name__ == "__main__":
    main('actual')