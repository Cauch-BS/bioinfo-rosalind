from dotenv import dotenv_values

def find_longest_path(
        all_edges:dict[int, dict[int, int]],
        source: int,
        max_val: int,
        sink: int
) -> tuple[int, list]:
    path_lengths: list[int] = [0 for _ in range(max_val + 1)]
    paths_memo: list[list] = [[] for _ in range(max_val + 1)]
    for node in all_edges:
        for vertex in all_edges[node]:
            if path_lengths[vertex] > path_lengths[node] + all_edges[node][vertex]:
                continue
            if source not in paths_memo[node] and node != source:
                continue
            paths_memo[vertex] = paths_memo[node] + [node]
            path_lengths[vertex] = path_lengths[node] + all_edges[node][vertex]
    paths_memo[sink].append(sink)
    return path_lengths[sink], paths_memo[sink]
            
def main():
    # Load the DNA sequence from the file
    datadir = dotenv_values(".env")["DATA_DIR"]
    loadfile = "rosalind_ba5d.txt"
    nucseq = open(f"{datadir}{loadfile}", "r+", encoding=None)

    # with nucseq as seq_file:
    #     _, _ = int(seq_file.readline()), int(seq_file.readline())
    #     lines = seq_file.readlines()
    #     seq_file.seek(2)
    #     for line in lines:
    #         line = line.strip()
    #         int1, int2, int3 = line.split()
    #         new_line = f"{int1}->{int2}:{int3}\n"
    #         seq_file.write(new_line)

    with nucseq as seq_file:
        source, sink = int(seq_file.readline()), int(seq_file.readline())
        max_val = sink
        all_edges = {node: dict() for node in range(sink)}
        for line in seq_file:
            if not line.strip():
                break
            edge, weight = line.strip().split(":")
            start, end = edge.split("->")
            start, end, weight = int(start), int(end), int(weight)
            if end > max_val:
                max_val = end
            try:
                all_edges[start][end] = weight
            except KeyError:
                all_edges[start] = {}
                all_edges[start][end] = weight

    result = find_longest_path(
        all_edges,
        source, 
        max_val,
        sink
    )

    print(result[0])
    print(*result[1], sep="->")

    return result

if __name__ == "__main__":
    main()
