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
basic_names = ["HC_", "SCH_"]
prefix = "example_"
new_names = "converted_" .* basic_names

# ===-===-===-===-===-
total_matrices = 44
all_matrices = Dict()
inversed_matrices = Dict()
normed_matrices = Dict()
data_keys = ["hc", "sch"]
for k in data_keys
    all_matrices[k] = Any[]
end

for f in 0:total_matrices-1, k in data_keys
    f_name = prefix * uppercase(string(k)) * "_"
    full_path = rawdatadir("original_weight_matrices", f_name * "$(f)" * ".csv")
    push!(all_matrices[k], readdlm(full_path, ',', Float64, '\n'))
end

for k in data_keys
    inversed_matrices[k] = inverse_elements(all_matrices[k])
end

# SAnity check- smallest elemnt in the upper diagonal should be largest in inver and vice versa
upper_diagonal_indices = [x[1] != x[2] for x in [k for k in CartesianIndices(all_matrices["hc"][1])]]
all_upper = [k for k in all_matrices["hc"][1][upper_diagonal_indices] if k != 0]
inversed_upper = [k for k in inversed_matrices["hc"][1][upper_diagonal_indices] if k != 0]
map(x -> x(all_upper), [findmin, findmax])
map(x -> x(inversed_upper), [findmin, findmax])

# ===-===-===-===-===-
for k in data_keys
    normed_matrices[k] = normalise_matrix(inversed_matrices[k])
end

# ===-===-===-===-===-
# Save matrices
prefix1 = "inversed_"
prefix2 = "inversed_normalised_"
extensions = ".csv"

ispath(prodatadir("inversed_matrices")) || mkpath(prodatadir("inversed_matrices"))

for k in data_keys
    for (index, mat) = enumerate(inversed_matrices[k])
        file_name = prefix1 * uppercase(string(k)) * "_$(index)" * extensions
        open(prodatadir("inversed_matrices", file_name), "w") do io
            writedlm(io, mat, ",")
        end
    end
end

for k in data_keys
    for (index, mat) = enumerate(normed_matrices[k])
        file_name = prefix2 * uppercase(string(k)) * "_$(index)" * extensions
        open(prodatadir("inversed_matrices", file_name), "w") do io
            writedlm(io, mat, ",")
        end
    end
end

