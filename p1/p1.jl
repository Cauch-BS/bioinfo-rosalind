using DotEnv
DotEnv.config()
datadir = ENV["DATA_DIR"]
open(joinpath(datadir, "rosalind_dna.txt")) do seqfile
    seq = readline(seqfile)
    result = []
    for base in ['A', 'C', 'G', 'T']
        push!(result, count(x -> x == base, seq), " ")
    end
    println(join(result))
end