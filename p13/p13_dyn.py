import numpy as np
from dotenv import dotenv_values

def find_minimal(
        coin_vec: np.ndarray,
        equal_to: int,
    ):
    memorize_mins = np.full((equal_to + 1), np.iinfo(np.int64).max)
    memorize_mins[0] = 0
    for m in range(1, equal_to + 1):
        for coin in coin_vec:
            if m >= coin:
                memorize_mins[m] = np.min([
                    memorize_mins[m],
                    memorize_mins[m - coin] + 1
                ])

    return memorize_mins[equal_to]

    
def main():
    # Load the DNA sequence from the file
    datadir = dotenv_values(".env")["DATA_DIR"]
    loadfile = "rosalind_ba5a.txt"
    data = open(f"{datadir}{loadfile}", "r", encoding=None)

    with data as data_file:
        equal_to = int(data_file.readline())
        equality_vec = np.array(list(
            map(int, data_file.readline().split(','))
        ), dtype = np.int64)


        print(
            find_minimal(
                coin_vec = equality_vec,
                equal_to = equal_to
            )
        )

if __name__ == "__main__":
    main()
