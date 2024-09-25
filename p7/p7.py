from itertools import product

from dotenv import dotenv_values
import numpy as np
from numpy.lib.stride_tricks import as_strided


def get_hamming_distance(
    strided_seq: np.ndarray, motif_candidates: np.ndarray
) -> np.ndarray:
    all_hamming = np.sum(strided_seq[:, np.newaxis] != motif_candidates, axis=2)
    min_hamming = np.min(all_hamming, axis=0)
    # print(f">>> Min Hamming Distances are {min_hamming}")
    return min_hamming


def main(mode: str = "debug"):
    # Load the DNA sequence from the file
    datadir = dotenv_values(".env")["DATA_DIR"]
    loadfile = "rosalind_ba2b.txt"
    nucseq = open(f"{datadir}{loadfile}", "r", encoding=None)
    nuc_to_inc = bytes.maketrans(b"ACGT", b"\x00\x01\x02\x03")
    nucleotides = np.array(list("ACGT"), dtype="U1")

    with nucseq as seq_file:
        motif_length = int(seq_file.readline().strip())
        strided_seqs = list()
        if "debug" in mode:
            seq_list = seq_file.readline().strip().split()
        elif "actual" in mode:
            seq_list = []
            for line in seq_file:
                seq_list.append(line.strip())
        else:
            raise ValueError(f"Invalid mode: {mode}")

        for line in seq_list:
            seq = np.frombuffer(
                line.strip().encode().translate(nuc_to_inc), dtype=np.uint8
            )
            seqlen = seq.shape[0]
            newshape = (seqlen - motif_length + 1, motif_length)
            newstrides = (seq.strides[0], seq.strides[0])
            strided_seq = as_strided(seq, shape=newshape, strides=newstrides)
            strided_seqs.append(strided_seq)

        motif_candidates = np.array(list(product(range(4), repeat=motif_length)))
        sum_hamming_distance = np.zeros(motif_candidates.shape[0])

        for strided_seq in strided_seqs:
            sum_hamming_distance += get_hamming_distance(
                strided_seq,
                motif_candidates,
            )
        motifs = nucleotides[motif_candidates[
            sum_hamming_distance == np.min(sum_hamming_distance)
        ]]
        motif_list = ["".join(row) for row in motifs]
        print(*motif_list)

if __name__ == "__main__":
    main("actual")
