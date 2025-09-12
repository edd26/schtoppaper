#=
Loads original weight matrix and converts the values such that the zero represents
no connection, small values represent weak connection, values closer to 1 represent
strong connection. Also, saves the matrices in the data folder.
=#
using DrWatson
@quickactivate "schtoppaper"
# ===-===-===-===-===-
using DelimitedFiles
using LinearAlgebra
using Statistics

using Plots

"MatrixUtils.jl" |> srcdir |> include
# ===-===-===-===-===-
rawdatadir(args...) = datadir("exp_raw", args...)
prodatadir(args...) = datadir("exp_pro", args...)

# ===-===-===-===-===-
basic_names = ["HCP_"]
prefix = "fdt_network_matrix_norm_"
new_names = "converted_" .* "fdt_network"

# ===-===-===-===-===-
total_matrices = 166
all_matrices = Dict()
inversed_matrices = Dict()
double_inversed_matrices = Dict()
normed_matrices = Dict()
data_keys = ["HCP"]
for k in data_keys
    all_matrices[k] = Any[]
end

## ===-===-===-===-===-===-
# Load matrices >>>
using MAT

raw_connectivity_matrices = Dict{String,Matrix{Float64}}()
all_subjects = rawdatadir("HCP", "raw_data") |> readdir

subject_folder = all_subjects[1]
for subject_folder in all_subjects
    try
        mat_file = rawdatadir("HCP", "raw_data", subject_folder, "1_AAL", "fdt_network_matrix.mat") |> matread
        raw_connectivity_matrices[subject_folder] = mat_file["SC"]
    catch
        @warn "File not found for subject: $(subject_folder)"
        continue
    end
end

all_matrices["HCP"] = [mat for (sub, mat) in raw_connectivity_matrices] |> vcat

# Load matrices <<<
## ===-===-===-===-===-===-
matrix_set = all_matrices["HCP"]
inversed_matrices["HCP"] = inverse_elements(all_matrices["HCP"])
double_inversed_matrices["HCP"] = inverse_elements(inversed_matrices["HCP"])

# SAnity check- smallest elemnt in the upper diagonal should be largest in inver and vice versa
off_diagonal_indices = [x[1] != x[2] for x in [k for k in CartesianIndices(all_matrices["HCP"][1])]]
all_upper = [k for k in all_matrices["HCP"][1][off_diagonal_indices] if k != 0]
inversed_upper = [k for k in inversed_matrices["HCP"][1][off_diagonal_indices] if k != 0]
map(x -> x(all_upper), [findmin, findmax])
map(x -> x(inversed_upper), [findmin, findmax])

findall(x -> x == 1, all_matrices["HCP"][1])

all_matrices["HCP"][1][27, 84]
all_matrices["HCP"][1][84, 27,]

# First question- is the data in the same format? the range of values for matrices in COBRE is (1e-10,1e-3), while for the HCP is <1e0, 1e8)

# ===-===-===-===-===-
for k in data_keys
    normed_matrices[k] = normalise_matrix(inversed_matrices[k])
end

# ===-===-===-===-===-
# Save matrices
new_prefix = "inversed_"
extensions = ".csv"
k = "HCP"
for k in data_keys
    prodatadir("inversed_matrices", k,) |> isdir || prodatadir("inversed_matrices", k) |> mkdir

    for (index, mat) = enumerate(normed_matrices[k])
        file_name = new_prefix * uppercase(string(k)) * "_$(index)" * extensions
        open(prodatadir("inversed_matrices", k, file_name), "w") do io
            writedlm(io, mat, ",")
        end
    end
end

new_prefix = "original_"
extensions = ".csv"
k = "HCP"
for k in data_keys
    prodatadir("original_matrices", k,) |> isdir || prodatadir("original_matrices", k) |> mkpath

    for (index, mat) = enumerate(double_inversed_matrices[k])
        file_name = new_prefix * uppercase(string(k)) * "_$(index)" * extensions
        open(prodatadir("original_matrices", k, file_name), "w") do io
            writedlm(io, mat, ",")
        end
    end
end
