from dotenv import dotenv_values
import numpy as np

def find_best(  # type: ignore
        alignment: np.ndarray,
        seq_1: np.ndarray,
        seq_2: np.ndarray,
        indel: int
    ) -> tuple[int, np.ndarray, np.ndarray]:
    n, m = alignment.shape
    print(f"Requires {n} * {m} iterations: {n * m}")
    n_1, m_1 = n + 1, m + 1
    memo_path_length = np.zeros((n_1, m_1), dtype=np.int32)
    memo_path = np.full((n_1, m_1), '', dtype="U1")

    for i in range(n_1):
        for j in range(m_1):
            if i == 0 or j == 0:
                continue  
            down_val = memo_path_length[i-1][j] - indel
            right_val = memo_path_length[i][j-1] - indel
            oblique_val = memo_path_length[i-1][j-1] + alignment[i-1][j-1]
            
            max_val = max(0, down_val, right_val, oblique_val) 
            
            if max_val == down_val:
                memo_path[i][j] = 'D'
            elif max_val == right_val:
                memo_path[i][j] = 'R'
            elif max_val == oblique_val:
                memo_path[i][j] = 'O'

            memo_path_length[i][j] = max_val

    i, j = np.unravel_index(np.argmax(memo_path_length), memo_path_length.shape)
    aligned_1 = []
    aligned_2 = []
    while i > 0 and j > 0 and memo_path_length[i][j] > 0:
        go_to = memo_path[i][j]
        if go_to == "D":
            i -= 1
            aligned_1.append(seq_1[i])
            aligned_2.append(-1)  
        elif go_to == "R":
            j -= 1
            aligned_1.append(-1)
            aligned_2.append(seq_2[j])
        elif go_to == "O":
            i -= 1
            j -= 1
            aligned_1.append(seq_1[i])
            aligned_2.append(seq_2[j])
    aligned_1.reverse()
    aligned_2.reverse()

    return int(np.max(memo_path_length)), np.array(aligned_1, dtype=np.int32), np.array(aligned_2, dtype=np.int32)

def main() -> None:
    PAM = np.array(
        [
            [  2, -2,  0,  0, -3,  1, -1, -1, -1, -2, -1,  0,  1,  0, -2,  1,  1,  0, -6, -3],
            [ -2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4,  0, -2, -2, -8,  0],
            [  0, -5,  4,  3, -6,  1,  1, -2,  0, -4, -3,  2, -1,  2, -1,  0,  0, -2, -7, -4],
            [  0, -5,  3,  4, -5,  0,  1, -2,  0, -3, -2,  1, -1,  2, -1,  0,  0, -2, -7, -4],
            [ -3, -4, -6, -5,  9, -5, -2,  1, -5,  2,  0, -3, -5, -5, -4, -3, -3, -1,  0,  7],
            [  1, -3,  1,  0, -5,  5, -2, -3, -2, -4, -3,  0,  0, -1, -3,  1,  0, -1, -7, -5],
            [ -1, -3,  1,  1, -2, -2,  6, -2,  0, -2, -2,  2,  0,  3,  2, -1, -1, -2, -3,  0],
            [ -1, -2, -2, -2,  1, -3, -2,  5, -2,  2,  2, -2, -2, -2, -2, -1,  0,  4, -5, -1],
            [ -1, -5,  0,  0, -5, -2,  0, -2,  5, -3,  0,  1, -1,  1,  3,  0,  0, -2, -3, -4],
            [ -2, -6, -4, -3,  2, -4, -2,  2, -3,  6,  4, -3, -3, -2, -3, -3, -2,  2, -2, -1],
            [ -1, -5, -3, -2,  0, -3, -2,  2,  0,  4,  6, -2, -2, -1,  0, -2, -1,  2, -4, -2],
            [  0, -4,  2,  1, -3,  0,  2, -2,  1, -3, -2,  2,  0,  1,  0,  1,  0, -2, -4, -2],
            [  1, -3, -1, -1, -5,  0,  0, -2, -1, -3, -2,  0,  6,  0,  0,  1,  0, -1, -6, -5],
            [  0, -5,  2,  2, -5, -1,  3, -2,  1, -2, -1,  1,  0,  4,  1, -1, -1, -2, -5, -4],
            [ -2, -4, -1, -1, -4, -3,  2, -2,  3, -3,  0,  0,  0,  1,  6,  0, -1, -2,  2, -4],
            [  1,  0,  0,  0, -3,  1, -1, -1,  0, -3, -2,  1,  1, -1,  0,  2,  1, -1, -2, -3],
            [  1, -2,  0,  0, -3,  0, -1,  0,  0, -2, -1,  0,  0, -1, -1,  1,  3,  0, -5, -3],
            [  0, -2, -2, -2, -1, -1, -2,  4, -2,  2,  2, -2, -1, -2, -2, -1,  0,  4, -6, -2],
            [ -6, -8, -7, -7,  0, -7, -3, -5, -3, -2, -4, -4, -6, -5,  2, -2, -5, -6, 17,  0],
            [ -3,  0, -4, -4,  7, -5,  0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2,  0, 10],
        ],
        dtype = np.int32
    )
    SIGMA = 5
    AA_TO_INC = bytes.maketrans(
        b"ACDEFGHIKLMNPQRSTVWY", 
        bytes(range(20))
    )
    AA = np.array(list("ACDEFGHIKLMNPQRSTVWY-"), dtype="U1")

    datadir = dotenv_values(".env")["DATA_DIR"]
    loadfile = "rosalind_ba5f.txt"
    outfile =  open("answers/rosalind_ba5f_p18.txt", "w", encoding = "UTF-8")
    aaseq = open(f"{datadir}{loadfile}", "r", encoding=None)
    with aaseq as seq_file:
        seq_1 = np.frombuffer(
            seq_file.readline().strip().encode().translate(AA_TO_INC), dtype=np.uint8
        )
        seq_2 = np.frombuffer(
            seq_file.readline().strip().encode().translate(AA_TO_INC), dtype=np.uint8
        )
    alignment = PAM[np.ix_(seq_1, seq_2)]
    results = find_best(
        alignment,
        seq_1,
        seq_2,
        indel = SIGMA
    )

    print(results[0], file = outfile)
    print(''.join(AA[results[1]]), file = outfile)
    print(''.join(AA[results[2]]), file = outfile)


if __name__ == "__main__":
    main()


