"""
    TopoDims(min_dim::Int, max_dim::Int)

Structure for minimal and maximal dimensions used in Eirene computations.

Asserts if `min_dim <= max_dim`.
"""
# @option "topo dims" struct TopoDims
struct TopoDims
    min_dim::Int
    max_dim::Int

    function TopoDims(min_dim::Int, max_dim::Int)
        if min_dim > max_dim
            throw(ArgumentError("Condition not met: min_dim<=max_dim"))
        end
        new(min_dim, max_dim)
    end
end


"""
    get_dim(topo_dims::TopoDims)

TopoDims general getter.
"""
function get_dim(topo_dims::TopoDims, which::Union{Symbol, String})
    if String(which) == "min"
        return topo_dims.min_dim
    elseif String(which) == "max"
        return topo_dims.max_dim
    else
        throw(ArgumentError, "Unrecognized second argument. Have to be `min` or `max`")
    end
end

get_min_dim(topo_dims::TopoDims) = get_dim(topo_dims, :min)
get_max_dim(topo_dims::TopoDims) = get_dim(topo_dims, :max)

