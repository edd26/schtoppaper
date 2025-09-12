using Pipe

include(srcdir("DataSource", "SymmetricMatrixSource.jl"))
include(srcdir("DataSource", "SignalMatrixLoader.jl"))
include("TopoParams.jl")

"""
    DataInfo(...


Structure with all information about used dataset.
- identifier- tag used to refer to this data in scripts. It is used as a key in
    dictionaries.
- samples_count- how many samples are there for this data; this information is taken from
    the SymmetricMatrixSource, therfore if SymmetricMatrixSource ensures it by returning how many
    points can be loaded, it is actual number of samples;
- 'preprocessing_pipeline' and 'samples_limiter' can be given as named arguemnts
"""
# @option
struct DataInfo
    # Required:
    identifier::String
    sym_matrix_source::DataSource
    topo_params::TopoParams

    # Optional and with default values
    samples_limiter::Int
    preprocessing_pipeline::Vector{Function}
    ordering_kwargs::NamedTuple

    function DataInfo(identifier::String,
        sym_matrix_source::DataSource,
        topo_params::TopoParams;
        samples_limiter::Union{Symbol,Int}=:none,
        preprocessing_pipeline::Vector{Function}=Function[],
        ordering_kwargs::NamedTuple=(assign_same_values=true,)
    )

        if samples_limiter == :none
            @debug "DataInfo: using all samples."
            samples_limiter = get_samples_count(sym_matrix_source)
        elseif typeof(samples_limiter) == Symbol
            "Unknown symbol used for limiter. Possible is: none." |> AssertionError |> throw
        else
            run_limit_assertions(samples_limiter, 1, get_samples_count(sym_matrix_source))
        end

        new(identifier,
            sym_matrix_source,
            topo_params,
            samples_limiter |> Int,
            preprocessing_pipeline,
            ordering_kwargs,
        )
    end

    function DataInfo(identifier::String,
        sym_matrix_source::DataSource,
        topo_params::TopoParams,
        samples_limiter::Int,
        preprocessing_pipeline::Vector{Function},
        ordering_kwargs::NamedTuple,
    )
        run_limit_assertions(samples_limiter, 1, get_samples_count(sym_matrix_source))

        new(identifier,
            sym_matrix_source,
            topo_params,
            samples_limiter,
            preprocessing_pipeline,
            ordering_kwargs,
        )
    end
end


function run_limit_assertions(samples_limiter, min_val, max_val)
    samples_limiter > max_val && ("Limiter can not be larger than total matrices." |>
                                  AssertionError |>
                                  throw)
    samples_limiter < min_val && ("Sample limiter must be a natural number." |> AssertionError |> throw)
end

"""
DataInfo getters.
"""
get_identifier(data_info) = data_info.identifier
get_matrix_source(data_info::DataInfo) = data_info.sym_matrix_source
get_topo_params(data_info::DataInfo) = data_info.topo_params
get_min_dim(data_info::DataInfo) = get_min_dim(data_info.topo_params)
get_max_dim(data_info::DataInfo) = get_max_dim(data_info.topo_params)
get_preprocessing(data_info::DataInfo) = data_info.preprocessing_pipeline
get_ordering_kwargs(data_info::DataInfo) = data_info.ordering_kwargs


function get_vals_range(vals_range::String)
    if occursin("-", vals_range,)
        min_val, max_val = split(vals_range, "-")
        @info min_val max_val
    else
        min_val = "1"
        max_val = vals_range
    end
    min_val, max_val = map(x -> parse(Int, x), [min_val, max_val])

    return min_val, max_val
end

function get_samples_range(data_info::DataInfo)
    min_val, max_val = get_vals_range(data_info.samples_limiter)
    return min_val:max_val
end
function get_samples_limit(data_info::DataInfo; only_max::Bool=true)
    return data_info.samples_limiter 
end


get_matrix_identifier(data_info::DataInfo) = data_info |> get_matrix_source |> get_matrix_identifier

function get_full_identifier(data_info::DataInfo)
    id_pt1 = data_info |> get_matrix_source |> get_matrix_identifier
    id_pt2 = data_info |> get_topo_params |> get_matrix_identifier
    return Dict("id" => join([id_pt1, id_pt2], ""))
end


## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
# exported from script 6c
function get_common_dims_range(data1::Dict{String,DataInfo})
    dims_ranges = hcat([TopologyPreprocessing.get_dims_range(val, ret=:vector) for (key, val) in data1]...)

    return findmin(dims_ranges[1, :])[1], findmax(dims_ranges[2, :])[1]
end

function get_dims_range(data1::DataInfo; ret_type::Symbol=:tuple)
    if ret_type == :tuple
        return get_min_dim(data1), get_max_dim(data1)
    end
    if ret_type == :vector
        return [get_min_dim(data1), get_max_dim(data1)]
    end
    if ret_type == :range
        return get_min_dim(data1):get_max_dim(data1)
    end
    error("Unknown ret_type value.")
end

function get_matrices(data_info::DataInfo; kwargs...)
    get_matrices(get_matrix_source(data_info); kwargs...)
end
