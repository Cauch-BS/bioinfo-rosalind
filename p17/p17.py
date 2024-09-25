from dotenv import dotenv_values
import numpy as np

def find_best(
        alignment: np.ndarray,
        seq_1: np.ndarray,
        seq_2: np.ndarray,
        indel: int
    ) -> tuple[int, np.ndarray, np.ndarray]:
    n, m = alignment.shape
    print(f"Requires {n} * {m} iterations: {n * m}")
    n_1, m_1 = n + 1, m + 1
    memo_path_length = np.zeros((n_1, m_1), dtype = np.int32)
    memo_traversed = np.full((n_1, m_1), False, dtype = bool)
    memo_path = np.full((n_1, m_1), '', dtype="U1")

    def update_path(i, j) -> None:
        if i > 0:
            down_val = memo_path_length[i-1][j] - indel
            if down_val > memo_path_length[i][j] or not memo_traversed[i][j]:
                memo_path_length[i][j] = down_val
                memo_path[i][j] = 'D'  
                memo_traversed[i][j] = True
        if j > 0:
            right_val = memo_path_length[i][j-1] - indel
            if right_val > memo_path_length[i][j] or not memo_traversed[i][j]:
                memo_path_length[i][j] = right_val
                memo_path[i][j] = 'R'
                memo_traversed[i][j] = True
        if i > 0 and j > 0:
            oblique_val = memo_path_length[i - 1][j - 1] + alignment[i - 1][j - 1]
            if oblique_val > memo_path_length[i][j] or not memo_traversed[i][j]:
                memo_path_length[i][j] = oblique_val
                memo_path[i][j] = 'O'
                memo_traversed[i][j] = True

    for j in range(m_1):
        for i in range(n_1):
            update_path(i, j)

    i, j, aligned_1, aligned_2 = n, m, [], []
    while i > 0 or j > 0:
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

    return int(memo_path_length[n][m]), np.array(aligned_1, dtype = np.int32), np.array(aligned_2, dtype = np.int32)

def main() -> None:
    BLOSUM = np.array(
        [
            [ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2],
            [ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
            [-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3],
            [-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2],
            [-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3],
            [ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3],
            [-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2],
            [-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1],
            [-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2],
            [-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1],
            [-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1],
            [-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2],
            [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3],
            [-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1],
            [-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2],
            [ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2],
            [ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2],
            [ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1],
            [-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2],
            [-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7],
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
    loadfile = "rosalind_ba5e.txt"
    outfile =  open("answers/rosalind_ba5e_p17.txt", "w", encoding = "UTF-8")
    aaseq = open(f"{datadir}{loadfile}", "r", encoding=None)
    with aaseq as seq_file:
        seq_1 = np.frombuffer(
            seq_file.readline().strip().encode().translate(AA_TO_INC), dtype=np.uint8
        )
        seq_2 = np.frombuffer(
            seq_file.readline().strip().encode().translate(AA_TO_INC), dtype=np.uint8
        )
    alignment = BLOSUM[np.ix_(seq_1, seq_2)]
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


