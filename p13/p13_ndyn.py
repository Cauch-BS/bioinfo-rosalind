import numpy as np
from scipy.optimize import linprog  #type:ignore
from dotenv import dotenv_values

def find_minimal(
        equality_vec: np.ndarray,
        equal_to: int,
        coefficients: np.ndarray
    ):
    """
    Cheap Function for Minimizing c.T@x subject to the condition that A@x = M. 
    Here, A, x, and c are all vectors (more mathematically correctly members of
    the module Z ⊕ Z ⊕ ... ⊕  Z with Z the ring of integers)
    """
    res = linprog(
        coefficients,
        A_eq = equality_vec,
        b_eq = equal_to,
        integrality = 1 # set to integers
    )

    if res.success:
        return int(np.sum(res.x))
    else:
        raise ValueError("Unable to find minimal solution")
    
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
        coefficients = np.ones(equality_vec.shape)

        print(
            find_minimal(
                equality_vec = equality_vec.reshape(1, -1),
                equal_to = equal_to,
                coefficients = coefficients
            )
        )

if __name__ == "__main__":
    main()
