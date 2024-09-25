import numpy as np
from dotenv import dotenv_values
from pyinstrument import Profiler

def find_longest(
        down_matrix: np.ndarray,
        right_matrix: np.ndarray
    ) -> int:
    n, m_1 = down_matrix.shape
    n_1, m = right_matrix.shape
    assert n_1 - n == 1, "Right Matrix must be of shape (n + 1) * m"
    assert m_1 - m == 1, "Down Matrix must be of shape n * (m + 1)"
    memo_path_length = np.zeros((n_1, m_1))

    def update_path(i, j) -> None:
        if i == 0 and j == 0:
            pass
        elif i == 0 and j != 0:
            memo_path_length[i][j] = memo_path_length[i][j - 1] + right_matrix[i][j - 1]
        elif i != 0 and j == 0:
            memo_path_length[i][j] = memo_path_length[i - 1][j] + down_matrix[i -1][j]
        else:
            memo_path_length[i][j] = max(
                memo_path_length[i][j - 1] + right_matrix[i][j - 1],
                memo_path_length[i -1][j] + down_matrix[i - 1][j]
            )
        
    for i in range(n_1):
        for j in range(m_1):
            update_path(i, j)

    return int(memo_path_length[n][m])

    
def main(loadfile):
    # Load the DNA sequence from the file
    datadir = dotenv_values(".env")["DATA_DIR"]
    # loadfile = "rosalind_ba5b.txt"
    data = open(f"{datadir}{loadfile}", "r", encoding=None)

    with data as data_file:
        n, _ = map(int, data_file.readline().split(" "))
        down_matrix = np.array([
                list(map(int, data_file.readline().split(" ")))
                for _ in range(n)
            ],
            dtype = np.int64
        )
        _ = data_file.readline()
        right_matrix = np.array([
                list(map(int, data_file.readline().split(" ")))
                for _ in range (n + 1)
            ],
            dtype = np.int64
        )
        print(
            find_longest(
                down_matrix = down_matrix,
                right_matrix = right_matrix
            )
        )

if __name__ == "__main__":
    main("rosalind_test.txt")