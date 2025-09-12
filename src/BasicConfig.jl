#=

Definitions of all possible configurations options

=#
# using Revise

using DrWatson
@quickactivate "SchiTopology"

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
using Configurations
# using TOML
using DelimitedFiles
using LinearAlgebra

using Pipe

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
configdir(args...) = datadir("config_files", args...)

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
include("TopoParams.jl")
include("DataInfo.jl")


## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
