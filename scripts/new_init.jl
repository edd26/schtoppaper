#=
Run eirene for selected symmetric matrices.

Every data must be described with 'DataInfo' structure in order to be used here.
Once that is done, the desired data 'identifier' have to be added to
'data_keys' found in this file.

Output of this script:
- symmetric_matrices- a dictionary with symmetric_matrices matrices put into a
    vector
- ord_matrices- a dictionary with ordered matrices put into a vector
- all_C_data- a dictionary with Topological reseults produced with eirene
    formed into a vector
- all_all_topo_features- a dictionary with topological features generated from
    'all_C_data' with eirene

Symmetric matrices can be loaded (e.g. connectivity matrices) or generated
    (e.g. geometric matrices)

=#
##

# TODO is the data created here the same as in the previous version of the scrip?

import DrWatson: @quickactivate, srcdir
@quickactivate "schtoppaper"

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
using Eirene
using TopologyPreprocessing

import Base.Threads: @spawn, nthreads, @sync

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
include(srcdir("helper_functions.jl"))
include(srcdir("ConfigIndividuals.jl"))
"OptionStructures.jl" |> srcdir |> include

# include(srcdir("BasicConfig.jl"))
# srcdir("DataSource", "MatrixLoader.jl") |> include

# using Main.DataConfig: data_info


## ===-===-===-===-===-===-===-===-===-===-===-===-===-
"ArgsParsing.jl" |> srcdir |> include
parsed_args = SchiArgPar.parse_clustering_commandline()

## ===-===-===-===-===-===-===-===-===-===-
# Default configs
@info "Processing arguments: " parsed_args
data_set = parsed_args["data_set"]
data_config = parsed_args["data_config"]

MIN_DIM = parsed_args["min_dim"]
MAX_DIM = parsed_args["max_dim"]

## ===-===-===-===-===-
# Check required arguments
data_set |> isnothing && ErrorException("Unrecognised or unspecified data set.") |> throw
data_config |> isnothing && ErrorException("Unspecified data set option number.") |> throw

data_info = configure_individuals_from_args(data_set, data_config, MIN_DIM, MAX_DIM)

## ===-===-===-===-===-===-===-===-===-===-===-===-===-
data_keys = [get_identifier(dat) for (key, dat) in data_info]

if "hc" in data_keys && "sch" in data_keys
    # data_keys =
    hc_position = findfirst(x -> x == "hc", data_keys)
    sch_position = findfirst(x -> x == "sch", data_keys)
    if sch_position < hc_position
        @info "Swapping..."
        data_keys[hc_position] = "sch"
        data_keys[sch_position] = "hc"
    end
end
# data_keys = all_identifiers

min_dims = [data_info[k] |> get_min_dim for k in data_keys]
max_dims = [data_info[k] |> get_max_dim for k in data_keys]
MIN_DIM = max(min_dims...)
MAX_DIM = min(max_dims...)

if (min(min_dims...) != max(max_dims...)) || min(max_dims...) != max(max_dims...)
    @warn "The dimensions specified across data info are not matching, only the smallest will be used"
end

dim_range = MIN_DIM:MAX_DIM

force_symmat_loading = parsed_args["force_computations"]
force_ordering = parsed_args["force_ordering"]
force_topology = parsed_args["force_topology"]
force_topo_features = parsed_args["force_topo_features"]
force_computations = parsed_args["force_computations"]

# with_step = parsed_args["with_step"]
normalisation_type = NormalisationType(parsed_args["normalisation_type"])

# TODO get normalised barcodes was changed- it has to be tested on two cases, with and without the step
# TODO after tested, subject 38 has to be checked for the differences in the results
# TODO all other results have to be re-run and checked 
## ===-===-===-===-===-===-===-===-===-===-===-===-===-
# Load data in a new way

if !@isdefined symmetric_matrices
    symmetric_matrices = Dict{String,Vector{Array{Int,2}}}()
    symmetric_matrices = populate_dict(Dict(), [data_keys])

    local_key = data_keys[1]
    for (index, local_key) in enumerate(data_keys)
        @info "Getting symmetric matrices for key: $(local_key)"
        selected_data = data_info[local_key]

        # total_matrices = get_samples_limit(selected_data)
        matrices_path = datadir("exp_pro", "symmetric_matrix", local_key)

        # for k in 1:total_matrices
        # TODO change this so that data is craeted when parrameters are changed
        # if local_key == "rand" || occursin("geom", local_key) || occursin("HCP", local_key)
        # total_matrices = selected_data |> get_matrix_source |> get_samples_count
        total_matrices = selected_data |> get_samples_limit
        local container = Dict("data_info" => selected_data, "key" => local_key, "total_matrices" => total_matrices)
        # else
        #     container = Dict("data_info" => selected_data)
        # end

        mat, s = produce_or_load(
            matrices_path,
            container, # container
            prefix="symmetric_matrix_$(local_key)", # prefix for savename
            tag=false, #github tag
            force=force_symmat_loading,
        ) do container
            mat1 = container["data_info"] |> get_matrix_source |> get_matrices
            preproc_funcs = selected_data |> get_preprocessing
            if haskey(container, "total_matrices")
                total_matrices = container["total_matrices"]
            else
                total_matrices = container["data_info"] |> get_matrix_source |> get_samples_count
            end

            mat = Dict(
                "symmetric_matrix" => [apply_processing(mat1[k, :, :], preproc_funcs) for k in 1:total_matrices]
            )
        end
        symmetric_matrices[local_key] = mat["symmetric_matrix"]#  matrices_collection
    end
    # end
else
    @info "Using dataset form workspace"
end


## ===-===-===-===-===-===-===-===-
# Ordering

@info "Starting ordering matrices..."
if !@isdefined ordered_matrices
    ordered_matrices = populate_dict!(Dict(), [data_keys,])
    # Dict()

    @sync for (index, local_key) in enumerate(data_keys)
        # for (index, local_key) in enumerate(data_keys)
        @info "Ordering for key: $(local_key)"
        selected_data = data_info[local_key]
        # ordered_matrices[key] = Any[]

        # Collect all arguments in a dictionary
        ordering_args_dict = Dict(
            "id" => data_info[local_key] |> get_matrix_identifier,
            "kwargs" => data_info[local_key] |> get_ordering_kwargs
        )

        @spawn for (mat_index, mat) in enumerate(symmetric_matrices[local_key])
            # for (mat_index, mat) in enumerate(symmetric_matrices[local_key])
            @info "\t matrix:$(mat_index)"
            ord_matrices_path = datadir("exp_pro", "ordered_matrices", "individuals", local_key)
            ordering_args_dict["k"] = mat_index
            ordering_args_dict["matrix"] = mat
            if ordering_args_dict["kwargs"][1] == :skip
                ordering_args_dict["ordering_type"] = ordering_args_dict["kwargs"][1]
            else
                if normalisation_type == to_unique_values::NormalisationType
                    ordering_args_dict["ordering_type"] = :cont
                elseif normalisation_type == sports::NormalisationType
                    ordering_args_dict["ordering_type"] = :with_step
                elseif normalisation_type == max_elements::NormalisationType
                    ordering_args_dict["ordering_type"] = :with_step
                else
                    ErrorException("Unrecognised options") |> throw
                end
            end

            # ord_mat, _ = produce_or_load(datadir("ord_matrices", local_key), # path
            ord_mat, _ = produce_or_load(ord_matrices_path,
                ordering_args_dict, # container
                get_ordered_matrix, # function
                prefix="matrix_ordering_$(local_key)", # prefix for savename
                tag=false, #github tag
                force=force_ordering,
            )
            push!(ordered_matrices[local_key], ord_mat["ordered_mat"])
        end
    end

    if isempty(keys(ordered_matrices))
        throw(DimensionMismatch("Loaded empty set of ordered matrices. Remove file or force execution."))
    end
    @info "Finished ordering."
else
    @info "Using ordered_matrices form workspace"
end

do_HCP_trim = false
if (("HCP_test_orig" in data_keys) || ("dHCP_test_orig" in data_keys)) && do_HCP_trim
    vals_above_diagonal = ceil(Int, 94 * (94 - 1) / 2)
    ordering_limit = vals_above_diagonal - 1958
    # ordering_limit = 1500
    # ordering_limit = 400
    # items_to_zero/vals_above_diagonal


    local_key = "HCP_test_orig"
    # mat = ordered_matrices[local_key][1]
    total_matrices = length(ordered_matrices[local_key])

    for m in 1:total_matrices
        mat = copy(ordered_matrices[local_key][m])
        items_over_limit = findall(x -> x >= ordering_limit, mat)
        mat[items_over_limit] .= ordering_limit

        ordered_matrices[local_key][m] = mat
    end
    # # # visualize capped ordered matrix
    # index = 6
    # f = Figure()
    # ax1 = CairoMakie.Axis(f[1, 1])
    # CairoMakie.heatmap!(ax1, ordered_matrices[local_key][index])

    # ax2 = CairoMakie.Axis(f[1, 2])
    # hmap = CairoMakie.heatmap!(ax2, ordered_matrices["hc"][index])
    # CairoMakie.Colorbar(f[1, end+1],)
    # f
end
findmax(ordered_matrices[local_key][1])


## ===-===-===-===-===-===-===-===-===-
# temporary section- should be incorporated to pipeline
# Get ordered matrices and lower the resolution
# original_ordered_matrices = copy(ordered_matrices)
# ordered_matrices= copy(original_ordered_matrices )
# ordered_matrix = ordered_matrices["hc_lowres200"][1]
# #
# for (index, key) in enumerate(data_keys)
#     @info "$(key)"
#     for (mat_index, mat) in enumerate(ordered_matrices[key])
#         @info "\tmat_index $(mat_index)"
#         lower_resolution!(mat, factor=200)
#     end
# end

# findall(x->x==60, ordered_matrices["hc_lowres"][1])
## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
# Topological analysis
@info "Starting the topology computations:"

if !@isdefined all_C_data
    # all_C_data = Dict{String, Vector{Dict{String, Any}}}()
    all_C_data = Dict()

    @sync for (index, local_key) in enumerate(data_keys)
        # for (index, local_key) in enumerate(data_keys)
        @info "Topo data for key: $(local_key)"
        all_C_data[local_key] = Any[]

        # Convert all arguments into single dictionary
        topo_args_dict = Dict()# "id"=>mat_identifier["id"],
        topo_args_dict["local_info"] = data_info[local_key]
        min_dim, max_dim = get_dims_range(data_info[local_key])

        if normalisation_type == to_unique_values::NormalisationType
            nothing
        elseif normalisation_type == sports::NormalisationType
            topo_args_dict["ordering_type"] = :with_step
        elseif normalisation_type == max_elements::NormalisationType
            topo_args_dict["ordering_type"] = :with_step
        else
            ErrorException("Unrecognised options") |> throw
        end

        if MAX_DIM == 4
            topo_args_dict["ord_type"] = :with_step
        end

        if (("HCP_test_orig" == local_key) || (("dHCP_test_orig" == local_key))) && do_HCP_trim
            @info "Doing the trimmed version..."
            topo_args_dict["do_HCP_trim"] = do_HCP_trim
            topo_args_dict["ord_limit"] = ordering_limit
        end

        @spawn for (mat_index, mat) in enumerate(ordered_matrices[local_key])
            # for (mat_index, mat) in enumerate(ordered_matrices[local_key])
            @info "\t matrix:$(mat_index)"
            topo_args_dict["k"] = mat_index
            topo_args_dict["matrix"] = mat
            eirene_results_path = datadir("exp_pro", "eirene_results", "individuals", local_key)

            params_dict = topo_args_dict
            # local_C, _ = produce_or_load(datadir("eirene_results", key), # path
            local_C, _ = produce_or_load(eirene_results_path,
                topo_args_dict, # container
                get_eirene_results_for_data, # function
                prefix="C_$(local_key)_max_dim=$(max_dim)", # prefix for savename
                tag=false, #github tag
                force=force_topology,
            )
            if local_C |> isempty
                throw(ErrorException("Empty eirene result found"))
            end
            push!(all_C_data[local_key], local_C["C"])
        end
    end

    @info "Completed the topology computations."
else
    @info "Using all_C_data form workspace"
end

# C1 = all_C_data[data_keys[1]][2]
# C2 = all_C_data[data_keys[2]][2]
# plotbarcode_pjs(C1, dim=0:3)
# plotbarcode_pjs(C2, dim=0:3)

## ===-===-
# Get topological features
if !@isdefined all_topo_features
    all_topo_features = Dict()

    # @sync for (index, local_key) in enumerate(data_keys)
    for (index, local_key) in enumerate(data_keys)
        @info "Topo features for key: $(local_key)"

        all_topo_features[local_key] = Any[]

        topo_args_dict = Dict()# "id"=>mat_identifier["id"],
        topo_args_dict["local_info"] = data_info[local_key]
        # topo_args_dict["max_dim"] = MAX_DIM
        # topo_args_dict["min_dim"] = MIN_DIM

        max_dim = get_max_dim(data_info[local_key])
        min_dim = get_min_dim(data_info[local_key])
        topo_results_path = datadir("exp_pro", "topo_features", "individuals", local_key)

        # TODO check if every matrix is loaded separately
        # @spawn for (mat_index, C) in enumerate(all_C_data[local_key])
        for (mat_index, C) in enumerate(all_C_data[local_key])
            # for (mat_index, C) in enumerate(all_C_data[key])
            @info "\t matrix:$(mat_index)"
            topo_args_dict["k"] = mat_index
            topo_args_dict["C"] = C

            if normalisation_type == to_unique_values::NormalisationType
                nothing
            elseif normalisation_type == sports::NormalisationType
                topo_args_dict["ordering_type"] = :with_step
            elseif normalisation_type == max_elements::NormalisationType
                topo_args_dict["ordering_type"] = :with_step
            else
                ErrorException("Unrecognised options") |> throw
            end

            if (("HCP_test_orig" == local_key) || (("dHCP_test_orig" == local_key))) && do_HCP_trim
                topo_args_dict["do_HCP_trim"] = do_HCP_trim
                topo_args_dict["ord_limit"] = ordering_limit
            end

            topo_features, _ = produce_or_load(
                topo_results_path, # path
                topo_args_dict, # container
                # get_topo_features,# function
                prefix="topo_features_$(local_key)_max_dim=$(max_dim)_min_dim=$(min_dim)", # prefix for savename
                tag=false, #github tag
                force=force_topo_features,
            ) do topo_args_dict
                local_info = topo_args_dict["local_info"]
                k = topo_args_dict["k"]
                C = topo_args_dict["C"]

                topo_params = get_topo_params(local_info)
                min_dim, max_dim = get_min_dim(topo_params), get_max_dim(topo_params)
                features_list = get_features(topo_params)

                return Dict("topo_features" => get_topo_features(
                    C,
                    features_list,
                    min_dim=min_dim,
                    max_dim=max_dim
                )
                )
            end

            push!(all_topo_features[local_key], topo_features["topo_features"])
        end # index
    end # key

    @info "Completed the topology descriptors."
else
    @info "Using all_topo_features form workspace"
end
