using DrWatson
@quickactivate "schtoppaper"

using DelimitedFiles
using Eirene
using PersistenceLandscapes
using Makie
using CairoMakie
using Pipe
using TopologyPreprocessing
using LinearAlgebra: triu

"ClusterMergePlottingMakie.jl" |> srcdir |> include
# ===-===-
angles = range(0, stop=2π, length=8 + 1)
x_coords = sin.(angles)
y_coords = cos.(angles)

function append_scaffold(ax, scaffold_matrix, x_coords, y_coords, colours_palette; data_matrix::Matrix{T}, lw::Number=0, alpha::Float64=1.0,) where {T<:Number}
    if lw == 0
        lw = 10
    end

    connections = [[k[1], k[2]] for k in @pipe findall(x -> x > 0, scaffold_matrix) |> filter(x -> x[1] < x[2], _)]
    total_connections = size(connections, 1)
    for row in 1:total_connections
        src, target = connections[row]# , :]
        colour_to_use = data_matrix[src, target]

        x_vals = [x_coords[src], x_coords[target]]
        y_vals = [y_coords[src], y_coords[target]]

        CairoMakie.lines!(
            ax,
            x_vals,
            y_vals,
            color=colours_palette[colour_to_use],
            linewidth=lw,
            alpha=alpha
        )
    end # row
    f
end

function miss_elements_on_diagonal(selected_mat)
    data_mat = copy(selected_mat)
    total_rows, _ = size(data_mat)
    data_mat = data_mat |> Matrix{Union{Int,Missing}}

    for r in 1:total_rows
        data_mat[r, r] = missing
    end
    return data_mat
end

# ===-===-
data_path(args...) = datadir("exp_raw", "high_dim_matrix", args...)
matrices_names = [k for k in data_path() |> readdir if k[1] != '.']


# process early_born
early_born = readdlm("early.csv" |> data_path, ';');
replace!(early_born, "x1" => 14, "x2" => 15)
early_born = early_born |> Matrix{Int}

# process late_born
late_born = readdlm("late4.csv" |> data_path, ';');
# late_born = readdlm("late.csv" |> data_path, ';');
replace!(late_born, "x1" => 25, "x2" => 26);
late_born = late_born |> Matrix{Int}


# process long_lived
long_lived = readdlm("long.csv" |> data_path, ';')
replace!(long_lived, 8 => 27)
replace!(long_lived, "x1" => 26, "x2" => 27);
long_lived = long_lived |> Matrix{Int}

# ===-===-
do_double_cycle = false


early_born =
    if do_double_cycle
        [0 1 9 3 22 11 14 15
            1 0 2 10 23 16 17 18
            9 2 0 4 24 19 20 21
            3 10 4 0 25 26 27 28
            22 23 24 25 0 5 12 8
            11 16 19 26 5 0 6 13
            14 17 20 27 12 6 0 7
            15 18 21 28 8 13 7 0]
    else
        early_born
    end

late_born =
    if do_double_cycle
        [0 1 27 13 24 10 2 3
            1 0 5 28 23 8 4 7
            27 5 0 12 22 11 6 14
            13 28 12 0 21 20 19 18
            24 23 22 21 0 9 25 17
            10 8 11 20 9 0 15 26
            2 4 6 19 25 15 0 16
            3 7 14 18 17 26 16 0]
    else
        late_born
    end

long_lived =
    if do_double_cycle
        [0 1 27 3 18 9 10 11
            1 0 2 28 19 12 13 14
            27 2 0 4 20 15 16 17
            3 28 4 0 21 22 23 24
            18 19 20 21 0 5 25 8
            9 12 15 22 5 0 6 26
            10 13 16 23 25 6 0 7
            11 14 17 24 8 26 7 0]
    else
        long_lived
    end
# ===-===-
# get
data_matrices = [early_born, late_born, long_lived] .|> TopologyPreprocessing.get_ordered_matrix

do_inverse = false
@info "How about inverse filtration? Set to $(do_inverse)"

if do_inverse

    for (k, _) in enumerate(data_matrices)
        mat = data_matrices[k]
        max_val = max(mat...)

        data_matrices[k] = .-mat .+ max_val

        total_rows = size(mat, 1)
        diagonal_indices = [CartesianIndex(k, k) for k in 1:total_rows]
        data_matrices[k][diagonal_indices] .= 0
    end
end
## ===-===-===-===-
# Get the landscapes
barcode_dim1(args...; kwargs...) = barcode(args...; dim=1, kwargs...)
function sort_by_birth(bd_data; kwargs...)
    cycles_sorting = sortperm(bd_data[:, 2])
    return bd_data[cycles_sorting, :]
end

barcodes_sorting =
    barcodes = map(x -> x |>
                        Matrix{Int} |>
                        (y -> eirene(y, maxdim=2)) |>
                        barcode_dim1 |>
                        sort_by_birth,
        data_matrices
    )
selected_cycle = 1
refactored_barcodes = [[k[selected_cycle, 1], k[selected_cycle, 2]] for k in barcodes]

# ===-===-
struct Thresholds
    low::Int
    high::Int
end

scaffolds_thresholds = [Thresholds(refactored_barcodes[k]...) for k in 1:3]

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
do_nothing = "ok"
