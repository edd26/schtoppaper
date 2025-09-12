include("TopoDims.jl")

"""
    TopoParams(...

Structure for arguments for eirene and for extracting topological properties.

"""
# @option "topo dims"
struct TopoParams
    topo_dims::TopoDims
    features_list::Vector{String} # TODO should be enum with values limmited to barcodes, bettis, and normed versions
    coordinates::Array # TODO coordinates should not be hold in parameters for topology computations
    coordinates_labels::Array

    function TopoParams(; min_dim::Int=1,
        max_dim::Int=1,
        features_list::Vector{String}=[],
        coordinates=Matrix[],
        coordinates_labels=Matrix[])
        new(TopoDims(min_dim, max_dim), features_list, coordinates, coordinates_labels)
    end

    function TopoParams(topo_dims::TopoDims,
        features_list::Vector{String},
        coordinates,
        coordinates_labels)
        new(topo_dims, features_list, coordinates, coordinates_labels)
    end

end

"""
    get_min_dim(topo_params::TopoParans)

Min dim getter.
"""
get_min_dim(topo_params::TopoParams) = get_dim(topo_params.topo_dims, :min)

"""
    get_max_dim(topo_params::TopoParans)

Max dim getter.
"""
get_max_dim(topo_params::TopoParams) = get_dim(topo_params.topo_dims, :max)

"""
    get_features(topo_params::TopoParams)

Featurs getter.
"""
get_features(topo_params::TopoParams) = topo_params.features_list



get_topo_dims(topo_params::TopoParams) = topo_params.topo_dims
get_coords(topo_params::TopoParams) = topo_params.coordinates
get_coords_labesls(topo_params::TopoParams) = topo_params.coordinates_labels

function get_matrix_identifier(topo_params::TopoParams)
    min_dim = get_min_dim(topo_params)
    max_dim = get_max_dim(topo_params)

    return "min_dim=$(min_dim)_max_dim=$(max_dim)"
end
