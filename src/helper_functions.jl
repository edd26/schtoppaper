using Base: nothing_sentinel
import DrWatson: @quickactivate, datadir, @unpack, srcdir
@quickactivate "schtoppaper"

"OptionStructures.jl" |> srcdir |> include

using Plots
using PlotThemes
using PlotUtils
using JSON
using DelimitedFiles
using Statistics
using Pipe
import LinearAlgebra: issymmetric, diag
using Distributions: Uniform

import TopologyPreprocessing.generate_random_matrix,
    TopologyPreprocessing.generate_random_point_cloud,
    TopologyPreprocessing.generate_geometric_matrix
# ===

resultdir(args...) = projectdir("data", "generated", args...)

# ===
"""
    populate_dict(input_dict, source_keys; final_structure=Any[])

A recursive function to create a dictionary within dictionaries.

Returns a sequence of dictionaries, in which keys are the elements of
`source_keys`.

"""
function populate_dict(input_dict, source_keys; final_structure=Any[])
    final_dictionary = copy(input_dict)
    dict_keys = copy(source_keys)
    populate_dict!(final_dictionary, dict_keys; final_structure=final_structure)

    return final_dictionary
end

function populate_dict!(final_dictionary, dict_keys; final_structure=Any[], dict_type=Dict)
    if length(dict_keys) == 0
        # Do final call
        final_dictionary = copy(final_structure)

        return final_dictionary
    else
        # do a recursive call with dict_keys reduced by 1
        first_keys = copy(dict_keys[1])
        other_keys = copy(dict_keys[2:end])

        a_dict = dict_type()
        for local_key in first_keys
            @debug "\t" local_key dict_keys[2:end]
            a_dict[local_key] = populate_dict!(final_dictionary,
                other_keys,
                final_structure=final_structure)
        end
        final_dictionary = a_dict
    end

    return final_dictionary
end

"""
    set_value_by_keys(plots_dictionary, all_keys, data)

A recursive function to set a value at given path of keys.

"""
function set_value_by_keys(plots_dictionary, all_keys, data)
    if do_copy
        local_dict = copy(plots_dictionary)
    else
        local_dict = plots_dictionary
    end

    if length(all_keys) == 1
        # Do final call
        local_dict[all_keys[1]] = data
    else
        # Do recursive call
        local_dict[all_keys[1]] = set_value_by_keys(local_dict[all_keys[1]], all_keys[2:end], data)
    end

    return local_dict
end

function set_value_by_keys!(local_dict, all_keys, data; do_copy=true)
    if length(all_keys) == 1
        # Do final call
        local_dict[all_keys[1]] = data
    else
        # Do recursive call
        local_dict[all_keys[1]] = set_value_by_keys!(local_dict[all_keys[1]], all_keys[2:end], data)
    end

    return local_dict
end
# ===

function symmetrize_matrix(matrix)
    new_matrix = copy(matrix)

    return (new_matrix' .+ new_matrix) ./ 2
end

function symmetrize_matrix_old(matrix)
    new_matrix = copy(matrix)
    total_rows, total_cols = size(new_matrix)
    for row = 1:total_rows, col = 1:total_cols
        new_val = new_matrix[row, col] + new_matrix[col, row]
        new_val /= 2
        new_matrix[row, col] = new_val
        new_matrix[col, row] = new_val
    end
    return new_matrix
end

function drop_half(matrix)
    new_matrix = copy(matrix)

    # Get values from above diagonal
    matrix_values = [UpperTriangular(new_matrix)...] |> unique
    total_values = length(matrix_values)

    # Get the middle value
    middle_of_value = total_values ÷ 3

    # Sortperm them
    # Drop the values below the middle value
    threshold_value = (matrix_values|>sort)[middle_of_value]
    new_matrix[new_matrix.>threshold_value] .= 0

    return new_matrix
end

"""

Applies Fisher transform to the values in the matrix, that is takes the
`atanh` of every element in the matrix.
"""
function z_transform_matrix(matrix::Matrix)
    new_matrix = atanh.(matrix)
    set_diagonal_to!(new_matrix, 0)
    new_matrix
end

function set_diagonal_to!(matrix::Matrix, number::Number)
    total_rows, total_cols = size(matrix)
    for (row, col) in zip(1:total_rows, 1:total_cols)
        matrix[row, col] = number
    end
    matrix
end

function reverse_sign(matrix::Matrix)
    return .-matrix
end

"""
Randomizes all the zero values in the matrix that are not on the diagonal
"""
function randomize_zeros(matrix::Matrix)
    zero_randomized_matrix = copy(matrix)

    zeros_positions_above_diagonal = [c for c in findall(x -> x == 0, zero_randomized_matrix) if c[1] > c[2]]

    if !isempty(zeros_positions_above_diagonal)
        max_in_matrix = findmax(matrix)[1]
        oper(args...) =
            if max_in_matrix < 0
                -(args...)
            else
                +(args...)
            end

        for c in zeros_positions_above_diagonal
            rand_number = rand(Uniform(0, 1))
            zero_randomized_matrix[c[1], c[2]] = oper(max_in_matrix, rand_number)
            zero_randomized_matrix[c[2], c[1]] = oper(max_in_matrix, rand_number)
        end
    end

    return zero_randomized_matrix
end

function prune_NaN(matrix::Matrix)
    findnans(args...) = findall(x -> isnan(x), args...)
    pruned_matrix = copy(matrix)

    nan_positions = [k[1] for k in matrix[1, :] |> findnans] |> unique
    for position in nan_positions
        pruned_matrix[:, position] .= 0
        pruned_matrix[position, :] .= 0
    end

    return pruned_matrix
end

"""
    set0_to_other_than(a_matrix::Matrix, sign::Symbol)

Return matrix which elements are in the range:
- (-∞, 0] if 'symbol' is set to ':-
- [0, ∞) if 'symbol' is set to ':+

Throws an 'ArgumentError' if 'symbol' is set to any other value.
"""
function set0_to_other_than(a_matrix::Matrix, sign::Symbol)
    matrix_duplicate = deepcopy(a_matrix)
    if sign == :+
        matrix_duplicate[matrix_duplicate.<0] .= 0
    elseif sign == :-
        matrix_duplicate[matrix_duplicate.>0] .= 0
    else
        ArgumentError("Unrecognised argument- 'sign' can only be ':+' or ':-'.") |> throw
    end

    return matrix_duplicate
end

leave_only_positive(a_matrix::Matrix) = set0_to_other_than(a_matrix, :+)
leave_only_negative(a_matrix::Matrix) = set0_to_other_than(a_matrix, :-)

# ===-===-===-===-===-

function generate_matrices(d::Dict)
    @unpack total_regions,
    cubeR_dims,
    total_matrices = d
    generated_structures_matrices = Dict()
    generated_structures_matrices["rand"] = [generate_random_matrix(total_regions)
                                             for x in 1:total_matrices]
    for cubeR_dim in cubeR_dims
        generated_structures_matrices["geom_R$(cubeR_dim)"] = [@pipe total_regions |>
                                                                     TopologyPreprocessing.generate_random_point_cloud(_, cubeR_dim) |>
                                                                     generate_geometric_matrix
                                                               for x in 1:total_matrices
        ]
    end

    return Dict("matrices" => generated_structures_matrices)
end

get_geometric_matrix(total_regions::Int, cubeR_dim::Int) = generate_random_point_cloud(total_regions, cubeR_dim) |>
                                                           generate_geometric_matrix
# ===


push_symmetric!(matrix_set; data_file=datadir()) =
    push!(matrix_set,
        symmetrize_matrix(readdlm(data_file, ',', Float64, '\n')))
# ===

function get_betti_plt_collection(ord_matrices, max_dim; data_name="some", verbouse=false, file_path="./results/", min_index=1)
    betti_collection = Any[]

    for k in min_index:10
        verbouse && @info "Currently working on " k
        C = eirene(ord_matrices[k], maxdim=max_dim, model="vr")
        bettis = get_bettis(C, max_dim)

        save(file_path * data_name * "_$(k)_max_dim$(max_dim).jld", "bettis", bettis)

        ref_betti = plot_bettis(bettis)
        push!(betti_collection, ref_betti)
    end
    return betti_collection
end

function get_betti_collection(ord_matrices, max_dim; data_name="some", verbouse=false, file_path="./results/", min_index=1)
    betti_collection = Any[]

    for k in min_index:10
        verbouse && @info "Currently working on " k
        C = eirene(ord_matrices[k], maxdim=max_dim, model="vr")
        bettis = get_bettis(C, max_dim)

        save(file_path * data_name * "_$(k)_max_dim$(max_dim).jld", "bettis", bettis)

        ref_betti = plot_bettis(bettis)
        push!(betti_collection, ref_betti)
    end
    return betti_collection
end



function get_regions_from_file(file_name::String)
    all_regions = JSON.parsefile(datadir("exp_raw", file_name)) # full path

    all_region_keys = keys(all_regions)
    all_regions_split = split.(all_region_keys)

    regions_pairs = Dict{Int,String}()

    for value in all_regions_split
        regions_pairs[parse(Int, value[1])] = value[2]
    end
    return regions_pairs
end

function get_region_coords(file_name::String)
    all_regions = JSON.parsefile(datadir("exp_raw", file_name)) # full path

    all_region_keys = keys(all_regions)
    all_regions_split = split.(all_region_keys)

    regions_coords = Dict{String,Vector{Float64}}()

    for value = all_regions_split
        related_vector = all_regions[join(value, " ")]
        regions_coords[value[2]] = related_vector
    end
    return regions_coords
end


function get_inversed_regions_from_file(file_name::String)
    all_regions = JSON.parsefile(datadir("exp_raw", file_name)) # full path

    all_region_keys = keys(all_regions)
    all_regions_split = split.(all_region_keys)

    regions_pairs = Dict{String,Int}()

    for value in all_regions_split
        regions_pairs[value[2]] = parse(Int, value[1])
    end
    return regions_pairs
end


"""
    load_schizophrenia_data(d::Dict)

Loads the default data set, after symmetrizing, with the parameters in d
(key - value):
- data_keys - vector with the kyes for resulting dictionary
- total_matrices - determines how many matrices to load
- prefix - prefix to data name

sample use:
all_matrices = load_schizophrenia_data(@dict all_my_keys total_matrices prefix)

"""
function load_schizophrenia_data(d::Dict)
    @unpack data_keys, total_matrices, prefix = d
    all_matrices = Dict{String,Vector{Array{Float64,2}}}()
    for k in data_keys
        all_matrices[k] = Any[]
    end

    # ===-===-
    if prefix == "inversed_"
        @info "Using rescaled data, with: 0 no connection, 1 strong connection."
    end # if
    for t in 1:total_matrices, k = data_keys
        data_file = datadir("exp_pro", "inversed_matrices", prefix * uppercase(string(k)) * "_$(t).csv")
        push_symmetric!(all_matrices[k], data_file=data_file)
    end

    return all_matrices
end

"""
    elementwise_average_matrix(matrices_set)

Creates a 3D matrix from the matrix set and then takes the average of all non-zero
values in that matrix.

"""
function elementwise_average_matrix(matrices_set)
    # convert matrix set to 3d tensor
    @debug size(matrices_set)
    total_matrices = size(matrices_set, 1)
    matrix_wid = size(matrices_set[1], 1)

    all_data = zeros(matrix_wid, matrix_wid, total_matrices)
    mean_matrix = zeros(matrix_wid, matrix_wid)
    std_matrix = zeros(matrix_wid, matrix_wid)

    for t in 1:total_matrices
        all_data[:, :, t] = matrices_set[t]
    end

    single_plane_indices = CartesianIndices(all_data[:, :, 1])
    diagonal_indices = findall(x -> x[1] == x[2], single_plane_indices)

    for ind in single_plane_indices
        all_values = all_data[ind[1], ind[2], :]
        non_zeros = findall(x -> x != 0, all_values)
        mean_matrix[ind] = mean(all_values[non_zeros])
        std_matrix[ind] = std(all_values[non_zeros])
    end

    mean_matrix[diagonal_indices] .= 0
    std_matrix[diagonal_indices] .= 0

    mean_matrix[findall(x -> isnan(x), mean_matrix)] .= 0
    std_matrix[findall(x -> isnan(x), std_matrix)] .= 0

    return mean_matrix, std_matrix
end


function get_ordered_matrix(params_dict::Dict)
    @unpack matrix, k, kwargs, id, ordering_type = params_dict
    @info ordering_type
    if ordering_type == :with_step
        return Dict("ordered_mat" => get_ordered_matrix2(matrix; with_step=true, kwargs...))
    elseif ordering_type == :skip

        return Dict("ordered_mat" => matrix)
    else
        return Dict("ordered_mat" => TopologyPreprocessing.get_ordered_matrix(matrix; kwargs...))
    end
end

function get_ordered_matrix2(in_matrix::Matrix;
    assign_same_values::Union{Bool,Symbol}=false,
    force_symmetry::Bool=false,
    small_dist_grouping::Bool=false,
    min_dist::Number=1e-16,
    total_dist_groups::Int=0,
    ordering_start::Int=1,
    with_step::Bool=true)


    # ==
    mat_size = size(in_matrix)
    ord_mat = zeros(Int, mat_size)

    # how many elements are above diagonal
    if issymmetric(in_matrix) || force_symmetry
        matrix_indices =
            generate_indices(mat_size, symmetry_order=true, include_diagonal=false)
        do_symmetry = true
    else
        matrix_indices = generate_indices(mat_size, symmetry_order=false)
        do_symmetry = false
    end
    total_elements = length(matrix_indices)

    # Collect vector of indices
    all_ind_collected = arr_to_vec(matrix_indices)

    # Sort indices vector according to inpu array
    index_sorting = sort_indices_by_values(in_matrix, all_ind_collected)

    ordering_number = ordering_start
    step_number = ordering_start
    for k = 1:total_elements
        next_sorted_pos = index_sorting[k]
        mat_ind = matrix_indices[next_sorted_pos]

        if k != 1
            prev_sorted_pos = index_sorting[k-1]
            prev_mat_ind = matrix_indices[prev_sorted_pos]
            cond1 = in_matrix[prev_mat_ind] == in_matrix[mat_ind]

            if assign_same_values

                cond2 = small_dist_grouping
                cond3 = abs(in_matrix[prev_mat_ind] - in_matrix[mat_ind]) < min_dist

                if cond1 || (cond2 && cond3)
                    ordering_number -= 1
                end
            end

            if with_step && !cond1
                ordering_number = step_number
            end
        end
        order = ordering_number

        set_values!(ord_mat, mat_ind, order; do_symmetry=do_symmetry)
        ordering_number += 1
        step_number += 1
    end

    if with_step
        max_val = findmax(ord_mat)[1]
        max_possible = total_elements

        if max_val != max_possible
            all_highest = findall(x -> x == max_val, ord_mat)
            ord_mat[all_highest] .= max_possible
        end
    end

    return ord_mat
end


function get_max_bd_values(barcodes_patient; selected_dim=1, keys_set=["hc" "sch"])
    m_birth = Dict()
    m_death = Dict()
    for key in keys_set
        m_birth[key] = findmax(barcodes_patient[key][selected_dim][:, 1])[1]
        m_death[key] = findmax(barcodes_patient[key][selected_dim][:, 2])[1]
    end

    # an automatic expression for undefined number of keys could be used here
    if m_birth["hc"] >= m_birth["sch"]
        max_birth = m_birth["hc"]
    else
        max_birth = m_birth["sch"]
    end

    if m_death["hc"] >= m_death["sch"]
        max_death = m_death["hc"]
    else
        max_death = m_death["sch"]
    end

    return max_birth, max_death
end

function get_eirene_results_for_data(params_dict::Dict)
    local_info = params_dict["local_info"]
    k = params_dict["k"]
    matrix = params_dict["matrix"]

    topo_params = get_topo_params(local_info)
    min_dim, max_dim = get_min_dim(topo_params), get_max_dim(topo_params)
    coordinates = get_coords(topo_params)
    coordinates_labels = get_coords_labesls(topo_params)

    if !isempty(coordinates_labels)
        C = eirene(matrix, maxdim=max_dim, model="vr", pointlabels=coordinates_labels)
    else
        C = eirene(matrix, maxdim=max_dim, model="vr")
    end

    if !isempty(coordinates)
        C["input"]["model"] = "pc"
        C["input"]["pc"] = "genera"
        C["input"]["genera"] = coordinates
    end

    return Dict("C" => C)
end

function local_get_normalised_barcodes(barcodes, betti_numbers, normalisation_type::NormalisationType)

    if normalisation_type == to_unique_values::NormalisationType
        final_total_steps = betti_numbers[:, 1] |> unique |> length
    elseif normalisation_type == sports::NormalisationType
        final_total_steps = betti_numbers[end-1, 1]
    elseif normalisation_type == max_elements::NormalisationType
        final_total_steps = betti_numbers[end, 1]
    else
        ErrorException("Unknown type of normalisation") |> throw
    end

    return barcodes ./ final_total_steps
end

function get_topo_features(C, features_list; min_dim::Int=1, max_dim::Int=3)
    @debug features_list # TODO feature list as an enum would be much better
    topo_features = Dict{String,Any}()

    if "norm_barcodes" in features_list
        @debug "\t\t\tDoing norm_barcodes"
        barcodes = get_barcodes(C, max_dim; min_dim=min_dim)
        betti_curve = Eirene.betticurve(C, dim=1)
        topo_features["norm_barcodes"] = local_get_normalised_barcodes(barcodes, betti_curve, NormalisationType(0))
    end

    if "sports_norm_barcodes" in features_list
        @debug "\t\t\tDoing sports_norm_barcodes"
        barcodes = get_barcodes(C, max_dim; min_dim=min_dim)
        betti_curve = Eirene.betticurve(C, dim=1)
        topo_features["sports_norm_barcodes"] = local_get_normalised_barcodes(barcodes, betti_curve, NormalisationType(1))
    end

    if "max_norm_barcodes" in features_list
        @debug "\t\t\tDoing max_norm_barcodes"
        barcodes = get_barcodes(C, max_dim; min_dim=min_dim)
        betti_curve = Eirene.betticurve(C, dim=1)
        topo_features["max_norm_barcodes"] = local_get_normalised_barcodes(barcodes, betti_curve, NormalisationType(2))
    end

    if "barcodes" in features_list
        @debug "\t\t\tDoing barcodes"
        topo_features["barcodes"] = get_barcodes(C, max_dim; min_dim=min_dim)
    end

    if "bettis" in features_list
        @debug "\t\t\tDoing bettis"
        topo_features["bettis"] = get_vectorized_bettis(C, max_dim; min_dim=min_dim)
    end

    return topo_features
end


function get_patient_topo_features(params_dict::Dict)
    """
    Dr watson wrapper for get_patient_topo_features.
    """
    if length(params_dict) == 5
        @info "At 5"
        @unpack C_data, t, key, max_dim, selected_dim = params_dict

        return get_patient_topo_features(C_data, max_dim, selected_dim)
    else
        @unpack C_data, t, key, min_dim, max_dim, selected_dim = params_dict
        @info "At 6"

        return get_patient_topo_features(C_data, max_dim, selected_dim; min_dim=min_dim)
    end
end

function get_patient_topo_features(C_data, max_dim, selected_dim; min_dim=1)
    topo_features = Dict{String,Any}()
    barcodes_patient = get_barcodes(C_data, max_dim; min_dim=min_dim)
    bettis_patient = get_vectorized_bettis(C_data, max_dim; min_dim=min_dim)

    topo_features["barcodes"] = barcodes_patient
    topo_features["bettis"] = bettis_patient

    topo_features["norm_barcodes"] =
        TopologyPreprocessing.get_normalised_barcodes(barcodes_patient,
            bettis_patient)

    return topo_features
end

function get_patient_second_degree_topo_features(params_dict::Dict)
    """
    Dr watson wrapper for get_patient_second_degree_topo_features.
    """
    @unpack topo_features, t, key, max_dim = params_dict

    return get_patient_second_degree_topo_features(topo_features)
end
function get_patient_second_degree_topo_features(topo_features::Dict{String,Any})
    b_area = TopologyPreprocessing.get_area_under_betti_curve(topo_features["bettis"])
    mx_bettis = TopologyPreprocessing.get_max_bettis(topo_features["bettis"])

    all_lifetimes = get_barcode_lifetime(topo_features["norm_barcodes"], max_dim=max_dim)
    mx_lifetime = get_barcode_max_lifetime(all_lifetimes)

    bd_ratios = get_birth_death_ratio(topo_features["norm_barcodes"], max_dim=max_dim)
    mx_bd_raio = get_barcode_max_db_ratios(bd_ratios)

    results_dict = Dict()
    results_dict["bettis_areas"] = b_area
    results_dict["bettis_max_values"] = mx_bettis
    results_dict["barcode_lifetimes"] = all_lifetimes
    results_dict["barcode_max_lifetimes"] = mx_lifetime
    results_dict["barcode_bd_ratio"] = bd_ratios
    results_dict["barcodes_max_db"] = mx_bd_raio

    return results_dict
end



function get_cardinality_distribution_plot(faces_patient_data, data_keys)
    total_matrices = length(faces_patient_data[data_keys[1]])
    cardinality_distribution = Dict()
    for key in data_keys
        cardinality_distribution[key] = zeros(total_matrices)
        for k = 1:total_matrices
            cardinality_distribution[key][k] = length(faces_patient_data[key][k])
        end
    end

    bins_range = 0:10:1.1*(findmax(cardinality_distribution[data_keys[1]])[1])
    plt_kwargs = (bins=bins_range, alpha=0.5)

    plt0 = histogram(cardinality_distribution[data_keys[1]]; label="hc", plt_kwargs...)
    histogram!(cardinality_distribution[data_keys[1]]; label="sch", plt_kwargs...)
    title!("Cardinality values distribution")
    ylabel!("Cardinality")

    return plt0, cardinality_distribution
end

function get_bar_size(C_patient, normed_barcodes, dim)
    total_bars = size(normed_barcodes[dim], 1)
    barsizes = Array{Int64}(undef, total_bars)
    for i = 1:total_bars
        barsizes[i] = length(C_patient["cyclerep"][dim+2][i])
    end
    return barsizes
end


## Mocking graph measures data

function get_graph_measures(; data_file="",
    graph_measures=5,
    total_regions=94,
    total_subjects=86,
    matrix_source="matlab"
)
    size_tuple = (graph_measures, total_regions, total_subjects)
    if data_file == ""
        graph_matrix = rand(graph_measures, total_regions, total_subjects)
    else
        @info "reading file"
        graph_matrix = readdlm(data_file, ',', Float64, '\n')

        if size(graph_matrix) != size_tuple
            if length(graph_matrix) != graph_measures * total_regions * total_subjects
                @error "Can not convert to " graph_measures total_regions total_subjects
            end

            if matrix_source == "matlab"
                @info "matlab file"
                total_rows = size(graph_matrix, 1)
                final_matrix = zeros(size_tuple)
                for row in 1:total_rows
                    final_matrix[row, :, :] = reshape(graph_matrix[row, :], (total_regions, total_subjects))
                end
                graph_matrix = final_matrix
            else
                graph_matrix = reshape(graph_matrix, size_tuple)
            end
        end
    end
    return graph_matrix
end


"""

Takes a symmetric, ordered matrix and reduces the total number of unique values
by a `factor` (default values is 2).

Designed to reduce the number of cliques build during TDA analysis.
"""
function lower_resolution(ordered_matrix; kwargs...)
    matrix_duplicate = deepcopy(ordered_matrix)
    lower_resolution!(matrix_duplicate; kwargs...)
    return matrix_duplicate
end

function lower_resolution!(ordered_matrix::Matrix{T} where {T}; factor=2)
    original_matrix = copy(ordered_matrix)
    all_indices = Any[]

    factorial_elements = 1:(factor-1)
    for k in factorial_elements
        target_indices = findall(x -> x % factor == k, original_matrix)
        all_indices = vcat(all_indices, target_indices,)
        ordered_matrix[target_indices] .+= 1
    end
end


"""
    apply_processing(data, funcs::Vector{Function})


Recursive call of  functions from `funcs` to a matrix `data`.

Until `funcs` is a vector with one functions, the function calls 
function at the first position in the vector (`func[1]`) on the
    recurscive call where first argument (`data`) is the same and 
    second argument is vector of functions starting with second position
    (`funcs[2:end]`). Final call is call of the only element from vector 
    on the first argument.


Throws `ErrorException` if the recursion condition is not met.
"""
function apply_processing(data::Matrix, funcs::Vector{Function})
    local_data = copy(data)

    if length(funcs) == 0
        return data
    elseif length(funcs) == 1
        final_func = funcs[1]
        final_val = final_func(local_data)
        return final_val
    else
        first_func = funcs[1]
        left_funcs = funcs[2:end]
        return apply_processing(local_data, left_funcs) |> first_func
    end
    ErrorException("Failed to met recursive conditions!") |> throw
end

function resolve_min_dim(min_dim, data_info, data_keys)
    data_info_min_dims = [data_info[local_key] |> get_min_dim for local_key in data_keys]
    data_info_minimal = min(data_info_min_dims...)
    if data_info_minimal != min_dim
        @warn "Minimal dimension set with arguments does not match minimal dimension in data loader: $(data_info_minimal )!=$(min_dim)"
    end
    final_min_dim = max(data_info_minimal, min_dim)
    return final_min_dim
end
