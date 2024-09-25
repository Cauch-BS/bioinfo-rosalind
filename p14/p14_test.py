import numpy as np
from dotenv import dotenv_values

def generate_random_dataset(n, m):
    down_matrix = np.random.randint(0, 10, size=(n, m + 1))
    right_matrix = np.random.randint(0, 10, size=(n + 1, m))
    return down_matrix, right_matrix

def save_to_file(filename, n, m, down_matrix, right_matrix):
    with open(filename, 'w') as f:
        f.write(f"{n} {m}\n")
        for row in down_matrix:
            f.write(" ".join(map(str, row)) + "\n")
        f.write("-\n")
        for row in right_matrix:
            f.write(" ".join(map(str, row)) + "\n")

def main():
    n = 1271  # Number of rows
    m = 4132  # Number of columns
    down_matrix, right_matrix = generate_random_dataset(n, m)
    datadir = dotenv_values(".env")["DATA_DIR"]
    testfile = "rosalind_test.txt"
    data = f"{datadir}{testfile}"
    save_to_file(data, n, m, down_matrix, right_matrix)

if __name__ == "__main__":
    main()