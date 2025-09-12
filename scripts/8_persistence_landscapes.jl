#=

Generate persistence landscapes using PersistenceLandscape toolbox.
Take the average persistence landscapes for both patietns groups.

=#
# using Revise
import DrWatson: @quickactivate, scriptsdir
@quickactivate "schtoppaper"

# ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
"new_init.jl" |> scriptsdir |> include


selected_barcodes = ""
if normalisation_type == to_unique_values::NormalisationType
    selected_barcodes = "norm_barcodes"
elseif normalisation_type == sports::NormalisationType
    selected_barcodes = "sports_norm_barcodes"
elseif normalisation_type == max_elements::NormalisationType
    selected_barcodes = "max_norm_barcodes"
else
    ErrorException("Unrecognised options") |> throw
end
# ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
using PersistenceLandscapes

# ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
filter_out_infinite = true

# Missing variables
MIN_DIM, MAX_DIM = data_info[data_keys[1]] |> get_dims_range
dim_range = MIN_DIM:MAX_DIM
total_matrices = all_topo_features[data_keys[1]] |> length

# ===-===-===-===-===-===-===-===-===-
all_landscapes = populate_dict!(Dict(),
    [
        data_keys,
        [t for t in MIN_DIM:MAX_DIM]
    ],
)

data_key = "sch"
selected_patient = 1
dim_index = 1
use_outer_layer_only = parsed_args["use_outer_layer_only"]

for data_key in data_keys, selected_patient = 1:total_matrices
    barcodes_patient_normed =
        all_topo_features[data_key][selected_patient][selected_barcodes]

    for (dim_index, selected_dim) = dim_range |> enumerate
        barcodes = barcodes_patient_normed[dim_index]

        if filter_out_infinite && selected_dim == 0
            @warn "Filtering out the infinite interval, dim index$(dim_index) and  t=$(selected_dim)"
            barcodes = barcodes[1:end-1, :]
        end

        bar = [MyPair(barcodes[k, 1], barcodes[k, 2]) for k in 1:size(barcodes, 1)]
        pair_barcodes = PersistenceBarcodes(bar, selected_dim)

        pl1 = PersistenceLandscape(pair_barcodes, selected_dim)

        plandscape = nothing
        if use_outer_layer_only
            plandscape = PersistenceLandscape([pl1.land[1]], pl1.dimension)
        else
            plandscape = pl1
        end
        push!(all_landscapes[data_key][selected_dim], plandscape)
    end
end

# ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
# Compute the average for both groups
landscapes_collection = populate_dict!(Dict(),
    [
        data_keys,
        [t for t in MIN_DIM:MAX_DIM]
    ],
)

average_landscape = populate_dict!(Dict(),
    [
        data_keys,
        [t for t in MIN_DIM:MAX_DIM]
    ],
)

std_landscapes = populate_dict!(Dict(),
    [
        data_keys,
        [t for t in MIN_DIM:MAX_DIM]
    ],
)

for data_key in data_keys, (dim_index, selected_dim) = MIN_DIM:MAX_DIM |> enumerate
    vector_for_avg = all_landscapes[data_key][selected_dim]
    vec_sapce_of_landscapes = VectorSpaceOfPersistenceLandscapes(vector_for_avg)
    average_land = PersistenceLandscapes.real_average(vec_sapce_of_landscapes)

    landscapes_collection[data_key][selected_dim] = vec_sapce_of_landscapes
    average_landscape[data_key][selected_dim] = average_land
    std_landscapes[data_key][selected_dim] = standardDeviation(landscapes_collection[data_key][selected_dim])
end
