using DelimitedFiles  # For reading the file
using LinearAlgebra    # For matrix operations
using DotEnv
DotEnv.config()

function read_data(file_path)
    open(file_path) do data_file
        n, m = parse.(Int, split(readline(data_file)))  # Read n and m
        down_matrix = [parse.(Int, split(readline(data_file))) for _ in 1:n]
        readline(data_file)  # Skip the separator line
        right_matrix = [parse.(Int, split(readline(data_file))) for _ in 1:(n + 1)]
        return n, m, hcat(down_matrix...), hcat(right_matrix...)
    end
end

function find_longest(down_matrix::Array{Int64, 2}, right_matrix::Array{Int64, 2})::Int
    m_1, n = size(down_matrix)
    m, n_1 = size(right_matrix)
    @assert n_1 - n == 1 "Right Matrix must be of shape (n + 1) * m"
    @assert m_1 - m == 1 "Down Matrix must be of shape n * (m + 1)"
    memo_path_length = zeros(Int64, m_1, n_1)

    function update_path(i::Int, j::Int)
        if i == 0 && j == 0
            return
        elseif i == 0 && j != 0
            memo_path_length[i + 1, j + 1] = memo_path_length[i + 1, j] + down_matrix[i + 1, j]
        elseif i != 0 && j == 0
            memo_path_length[i + 1, j + 1] = memo_path_length[i, j + 1] + right_matrix[i, j + 1]
        else
            memo_path_length[i + 1, j + 1] = max(
                memo_path_length[i + 1, j] + down_matrix[i + 1, j],
                memo_path_length[i, j + 1] + right_matrix[i, j + 1]
            )
        end
    end

    for i in 0:m
        for j in 0:n
            update_path(i, j)
        end
    end

    return memo_path_length[m_1, n_1]
end

function main()
    datadir = ENV["DATA_DIR"] 
    loadfile = "rosalind_test.txt"
    n, m, down_matrix, right_matrix = read_data(datadir * loadfile)
    result = find_longest(down_matrix, right_matrix)
    println(result)
end

main()
