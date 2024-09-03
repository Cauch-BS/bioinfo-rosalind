using DotEnv
DotEnv.config()
datadir = ENV["DATA_DIR"]
open(joinpath(datadir, "rosalind_rna.txt")) do seqfile
    seq = readline(seqfile)
    println(replace(seq, "T" => "U"))
end