import DrWatson: srcdir, @strdict, @unpack
@quickactivate "SchiTopology"
"DataSource.jl" |> srcdir |> include

"""
Abstract type to bind symmetric matrix sources.
"""
abstract type SymmetricMatrixSource <: DataSource end


