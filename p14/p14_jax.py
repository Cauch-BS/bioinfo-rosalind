import numpy as np
from dotenv import dotenv_values
import jax
import jax.numpy as jnp

def find_longest(
        down_matrix: np.ndarray,
        right_matrix: np.ndarray
    ) -> int:
    n, m_1 = down_matrix.shape
    n_1, m = right_matrix.shape
    assert n_1 - n == 1, "Right Matrix must be of shape (n + 1) * m"
    assert m_1 - m == 1, "Down Matrix must be of shape n * (m + 1)"
    memo_path_length = jnp.zeros((n_1, m_1))
    indices = jnp.array(
        [[i, j] for i in range(n_1) for j in range(m_1)]
    )
    @jax.jit
    def main_body_fn(right_matrix, down_matrix, memo_path_length, indices):
        def update_path(carry, indices):
            i, j = indices
            def update_right():
                return jax.lax.select(
                    j >= 1,
                    carry[i, j - 1] + right_matrix[i , j - 1],
                    0.0
                )
            def update_down():
                return jax.lax.select(
                    i >= 1,
                    carry[i -1, j] + down_matrix[i - 1, j],
                    0.0
                )
            carry = carry.at[i, j].set(jnp.maximum(
                update_right(),
                update_down(),
            ))
            return carry, carry
        memo_path_length, _ = jax.lax.scan(update_path, memo_path_length, indices)
        return memo_path_length

    return int(main_body_fn(
        right_matrix,
        down_matrix,
        memo_path_length,
        indices
    )[n][m])

    
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