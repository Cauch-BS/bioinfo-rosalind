from dotenv import dotenv_values
import numpy as np
from numpy.lib.stride_tricks import as_strided


def main():
    # Load the DNA sequence from the file
    datadir = dotenv_values(".env")["DATA_DIR"]
    loadfile = "rosalind_ba2c.txt"
    nucseq = open(f"{datadir}{loadfile}", "r", encoding=None)
    nuc_to_inc = bytes.maketrans(b"ACGT", b"\x00\x01\x02\x03")
    nucleotides = np.array(list("ACGT"), dtype="U1")

    with nucseq as seq_file:
        seq = np.frombuffer(
            seq_file.readline().strip().encode().translate(nuc_to_inc), dtype=np.uint8
        )
        k_mer_len = int(seq_file.readline())
        profile_matrix = np.array(
            [list(map(float, line.split(" "))) for line in seq_file.readlines()],
            dtype=np.float32,
        )
        seqlen = seq.shape[0]
        newshape = (seqlen - k_mer_len + 1, k_mer_len)
        newstrides = (seq.strides[0], seq.strides[0])
        strided_seq = as_strided(seq, shape=newshape, strides=newstrides)
        profiled_seq = profile_matrix[strided_seq, np.arange(strided_seq.shape[1])]
        print(
            "".join(nucleotides[strided_seq[np.argmax(np.prod(profiled_seq, axis=1))]])
        )


if __name__ == "__main__":
    main()
