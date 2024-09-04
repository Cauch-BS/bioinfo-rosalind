# using DotEnv
# DotEnv.config()
# datadir = ENV["DATA_DIR"]
# open(joinpath(datadir, "inputs", "input_1.txt")) do seqfile
#     # Create a mapping from nucleotides to numbers
#     nucleotides = (UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T'))
#     nuc_to_inc = x -> findfirst(==(x), nucleotides)
#     seqline = strip(readline(seqfile))
#     seq = nuc_to_inc.(collect(UInt8.(codeunits(seqline))))
#     windowline = strip(readline(seqfile))
#     window = nuc_to_inc.(collect(UInt8.(codeunits(windowline))))
#     seqlen = length(seq)
#     windowsize = length(window) - 1
#     stride_through_seq = reshape(
#         view(seq, 1 : seqlen - windowsize), 1, seqlen - windowsize
#     ) .+ (0 : windowsize)
#     println(stride_through_seq)
#     print(window)
# end