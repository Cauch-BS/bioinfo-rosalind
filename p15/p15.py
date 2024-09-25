from dotenv import dotenv_values


def main():
    # Load the DNA sequence from the file
    datadir = dotenv_values(".env")["DATA_DIR"]
    loadfile = "rosalind_ba5c.txt"
    nucseq = open(f"{datadir}{loadfile}", "r", encoding=None)

    with nucseq as seq_file:
        seq_1 = seq_file.readline().strip()
        seq_2 = seq_file.readline().strip()
        n, m = len(seq_1), len(seq_2)
        lcs_table = [
            [0 for _ in range(m + 1)]
            for _ in range(n + 1)
        ]
        for i in range(1, n +1):
            for j in range(1, m + 1):
                if seq_1[i - 1] == seq_2[j - 1]:
                    lcs_table[i][j] = lcs_table[i - 1][j - 1] + 1
                else:
                    lcs_table[i][j] = max(
                        lcs_table[i - 1][j],
                        lcs_table[i][j - 1]
                    )
        lcs = []
        while n > 0 and m > 0:
            if seq_1[n - 1] == seq_2[m - 1]:
                lcs.append(seq_1[n - 1])
                n -= 1
                m -= 1
            elif lcs_table[n - 1][m] > lcs_table[n][m - 1]:
                n -= 1
            else:
                m -= 1
        print(''.join(lcs[::-1]))

if __name__ == "__main__":
    main()
