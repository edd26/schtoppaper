using Clustering

## ===-===-===-
# Structures
"""
    ClusterSplit

Data structure designed to keep track of the branches splits for splitting the
cluster n-steps from the top of the cluster.
"""
struct ClusterSplit
    name::String
    range::UnitRange
    deph::Int
    height::Float64
end

function get_cluster_split_name(clust_split::ClusterSplit)
    return parse(Int, clust_split.name, base=2)
end

"""
    HclustExtended
    
Data structure `Hclust` from `Clustering` module extended with a vector
containing information on which tree leaves is the merge spanning.
"""
struct HclustExtended
    heights::Vector{Float64}
    labels::UnitRange{Int}
    linkage::Symbol
    merges::Matrix{Int}
    method::Symbol
    order::Vector{Int}
    ranges

    function HclustExtended(hclust::Hclust{T}) where {T<:Real}
        ranges = get_merges_ranges(hclust)

        new(
            hclust.heights,
            hclust.labels,
            hclust.linkage,
            hclust.merges,
            hclust.method,
            hclust.order,
            ranges)
    end
end

"""
    HclustMerged

Data structure to represent cluster merging with the n-times split from top.
"""
struct HclustMerged
    heights::Vector{Float64}
    labels::UnitRange{Int}
    linkage::Symbol
    merges::Matrix{Int}
    method::Symbol
    order::Vector{Int}
    ranges
    clusters
    merged_leaves

    function HclustMerged(hclust_e::HclustExtended, clusters, merged_leaves)
        new(
            hclust_e.heights,
            hclust_e.labels,
            hclust_e.linkage,
            hclust_e.merges,
            hclust_e.method,
            hclust_e.order,
            hclust_e.ranges,
            clusters,
            merged_leaves
        )
    end
end

## ===-===-===-
# Functions 
"""
    get_merges_ranges(selected_cluster::Hclust)::Vector{UnitRange{Int}}

Processes hierachical clustering data structure `selected_cluster.merges`
and returns a vector with range of leaves (in final cluster) that are
related to the `meges`.

"""
function get_merges_ranges(selected_cluster::Hclust)::Vector{UnitRange{Int}}

    total_rows = size(selected_cluster.merges, 1)
    leaves_range_per_merge = fill(0:0, (total_rows,))

    index = row = 478
    # 330 leaves only merges
    # 478 first double branch merge
    # for (index, row) = 1:477 |> enumerate

    for (index, row) = 1:total_rows |> enumerate
        leaf_val1, leaf_val2 = [0, 0]
        leaves_range1 = leaves_range2 = 0:0

        left_val, right_val = selected_cluster.merges[row, :]

        if left_val < 0
            original_data_position1 = abs(left_val)
            # translate from original data postion into leaves in the final cluster
            leaf_val1 = findall(x -> x == original_data_position1, selected_cluster.order)[1]
        end

        if right_val < 0
            original_data_position2 = abs(right_val)
            # translate from original data postion into leaves in the final cluster
            leaf_val2 = findall(x -> x == original_data_position2, selected_cluster.order)[1]
        end

        # if the values are positive, get the already existing range the branch are  covering
        if right_val > 0
            leaves_range1 = [k for k in leaves_range_per_merge[right_val]]
        end
        if left_val > 0
            leaves_range2 = [f for f in leaves_range_per_merge[left_val]]
        end

        # Sort out what are the ranges of values (0 is not used in merges structure)
        left_range_limit = findmin([k for k in [leaves_range1..., leaves_range2..., leaf_val1, leaf_val2] if k != 0])[1]
        right_range_limit = findmax([leaves_range1..., leaves_range2..., leaf_val1, leaf_val2])[1]

        # Get range
        leaves_range_per_merge[row] = left_range_limit:right_range_limit
    end
    return leaves_range_per_merge
end


"""
    find_superset_from_clusters(leaves_range::Vector{UnitRange{Int}}, clusters::Dict)

Finds the superset of the `leaves_range` in `clusters` and returns its name.
"""
function find_superset_from_clusters(leaves_range::UnitRange{Int}, clusters::Union{OrderedDict,Dict})
    cluster_ranges = [val.range for (key, val) in clusters]
    cluster_names = [key for (key, val) in clusters]

    cluster_position = findall(x -> issubset(leaves_range, x), cluster_ranges)

    if cluster_position |> isempty
        return ""
    else
        return cluster_names[cluster_position][1]
    end
end

"""
get_popularity(mat) = mat |>
                      sum(_, dims=2) |>
                      binarize_matrix |>
                      findall(x -> x != 0, _) |>
                      length

Compute popularity of a cycle by counting how many non-zero elements are there
in the sum of cycles occurance.
"""
get_popularity(mat) = @pipe mat |>
                            sum(_, dims=2) |>
                            binarize_matrix |>
                            findall(x -> x != 0, _) |>
                            length

"""
    split_cluster_n_times(leaves_range_per_merge, leaves_heights, target_deph)

Merges cluster branches `target_deph` times thus creating `2^target_deph`
clusters.

Each of the cluster merges are created at differnt height, but all of them are
`total_deph` branch-splitts from the top of the cluster (unless there it was
not possible to create that deph due to not enough splits in the branch).

"""
function split_cluster_n_times(
    extended_cluster,
    max_deph;
    min_height::Number=0,
    min_popularity::Int=0,
    ppl_presence::Matrix{T}=Array{Int}(undef, 0, 0),
    allow_final_splits::Bool=true,
    max_width=100,
    min_width=0
) where {T<:Number}
    leaves_heights = extended_cluster.heights[end:-1:1]
    inv_leaves_range_per_merge = extended_cluster.ranges[end:-1:1]

    clusters_merges = Dict("0" => ClusterSplit("0", inv_leaves_range_per_merge[1], 0, leaves_heights[1]))


    for (iteration, (child_leaves_range, cutting_height)) in zip(inv_leaves_range_per_merge[2:end], leaves_heights[2:end]) |> enumerate
        cluster_superset_name = find_superset_from_clusters(child_leaves_range, clusters_merges |> sort)

        parent_cluster = clusters_merges[cluster_superset_name]
        parent_width = length(parent_cluster.range)

        if parent_cluster.range == child_leaves_range
            continue
        else
            @debug "New cluster to process"
        end

        # This handles popularity of parent cluster
        if !isempty(ppl_presence)
            parent_cluster_popularity = (ppl_presence|>binarize_matrix)[:, parent_cluster.range] |>
                                        get_popularity
        else
            parent_cluster_popularity = Inf
        end

        # mid_boundary should always be left region oriented
        if child_leaves_range[1] == parent_cluster.range[1]
            left_child_bounds = child_leaves_range
            right_child_bounds = (child_leaves_range[end]+1):parent_cluster.range[end]
        elseif child_leaves_range[end] == parent_cluster.range[end]
            left_child_bounds = parent_cluster.range[1]:(child_leaves_range[1]-1)
            right_child_bounds = child_leaves_range
        else
            continue
            ErrorException("This case was not programmed. Contact the developer") |> throw
        end

        # ===-
        left_popularity = @pipe (ppl_presence|>binarize_matrix)[:, left_child_bounds] |>
                                get_popularity
        right_popularity = @pipe (ppl_presence|>binarize_matrix)[:, right_child_bounds] |>
                                 get_popularity

        # is the depth of cluster reached?
        skip_split = false

        if parent_cluster.deph >= max_deph
            @debug "Depth reached for cluster $(parent_cluster.name)"
            skip_split = true
        elseif cutting_height <= min_height
            @debug "Height criterium reached for cluster $(parent_cluster.name)"
            skip_split = true
        elseif min_popularity >= parent_cluster_popularity
            @info "Popularity criterium reached for cluster $(parent_cluster.name) with $(parent_cluster_popularity)- it won't be split further."
            skip_split = true
            if parent_width > max_width
                @info "\t but the width is too big"
                skip_split = false
            end


        elseif !allow_final_splits && @pipe ([left_popularity, right_popularity] .< min_popularity) |> Vector{Bool} |> any
            @debug "Popularity criterium reached if splitting $(parent_cluster.name) with $(left_popularity) and $(right_popularity)- it won't be split further."
            skip_split = true

        end

        if (@pipe (([left_child_bounds, right_child_bounds] .|> length) .> max_width) |> Vector{Bool} |> any) && !(parent_cluster.deph >= max_deph)
            skip_split = false
        end
        if (@pipe (([left_child_bounds, right_child_bounds] .|> length) .<= min_width) |> Vector{Bool} |> any) && !(parent_cluster.deph >= max_deph)
            skip_split = true
        end

        if !skip_split
            @debug "Deppening the cluster"

            # Create new clusters with new ranges
            clust = pop!(clusters_merges, parent_cluster.name)
            prefix = clust.name
            left_branch_name = prefix * "0"
            clusters_merges[left_branch_name] = ClusterSplit(left_branch_name, left_child_bounds, clust.deph + 1, cutting_height)
            right_branch_name = prefix * "1"
            clusters_merges[right_branch_name] = ClusterSplit(right_branch_name, right_child_bounds, clust.deph + 1, cutting_height)
        end # if splitting

        if all(i -> i == max_deph, [clust.deph for (key, clust) in clusters_merges])
            break
        else
            @debug "Not breaking early"
        end
    end

    # Make sure all names have the same length
    clusters_names = [k for k in clusters_merges |> keys]
    clusters_names_len = [k |> length for k in clusters_names]

    # is correction needed?
    if findmin(clusters_names_len)[1] != findmax(clusters_names_len)[1]
        target_length = findmax(clusters_names_len)[1]

        for k in clusters_names
            name_len = k |> length
            if name_len != target_length
                clust = pop!(clusters_merges, k)
                fixed_name = clust.name
                for n in 1:(target_length-name_len)
                    fixed_name *= "0"
                end
                clusters_merges[fixed_name] = ClusterSplit(fixed_name, clust.range, clust.deph, clust.height)
            end
        end
    end # if different

    clusters_merges = clusters_merges |> sort
    return clusters_merges
end

function translate_to_merges(selected_cluster, min_clusters, cluster_height)
    merged_cluster = cutree(selected_cluster, k=min_clusters, h=cluster_height)

    unique_names = merged_cluster |> unique |> sort
    clusters_merges = OrderedDict()

    longest_name = [length(string(k, base=2)) for k in unique_names] |> maximum

    for k in unique_names
        clust_name = string(k, base=2)
        clust_name = repeat("0", longest_name - length(clust_name)) * clust_name


        vals_vec = findall(x -> x == k, merged_cluster |> sort)
        clust_range = minimum(vals_vec):maximum(vals_vec)

        clusters_merges[clust_name] = ClusterSplit(clust_name, clust_range, 0, cluster_height)

    end

    # Make sure all names have the same length
    clusters_names = [k for k in clusters_merges |> keys]
    clusters_names_len = [k |> length for k in clusters_names]

    # is correction needed?
    if findmin(clusters_names_len)[1] != findmax(clusters_names_len)[1]
        target_length = findmax(clusters_names_len)[1]

        for k in clusters_names
            name_len = k |> length
            if name_len != target_length
                clust = pop!(clusters_merges, k)
                fixed_name = clust.name
                for n in 1:(target_length-name_len)
                    fixed_name *= "0"
                end
                clusters_merges[fixed_name] = ClusterSplit(fixed_name, clust.range, clust.deph, clust.height)
            end
        end
    end # if different

    clusters_merges = clusters_merges |> sort
    return clusters_merges
end


function translate_to_merges2(selected_cluster, min_clusters, cluster_height; do_name_replace::Bool=false)
    merged_cluster = cutree(selected_cluster, k=min_clusters, h=cluster_height)
    # Where are clusters?
    merged_cluster2 = merged_cluster[selected_cluster.order]

    unique_names = merged_cluster2 |> unique
    sorted_unique_names = unique_names |> sort

    longest_name = [length(string(k, base=2)) for k in unique_names] |> maximum

    clusters_merges = OrderedDict()
    for k in sorted_unique_names
        clust_name = string(k, base=2)
        clust_name = repeat("0", longest_name - length(clust_name)) * clust_name

        vals_vec = findall(x -> x == k, merged_cluster2)
        clust_range = minimum(vals_vec):maximum(vals_vec)

        clusters_merges[clust_name] = ClusterSplit(clust_name, clust_range, 0, cluster_height)
    end
    # Make sure all names have the same length
    clusters_names = [k for k in clusters_merges |> keys]
    clusters_names_len = [k |> length for k in clusters_names]

    # ===-
    # Replace names in the dictionary and in the clust split
    # replace names if sorted are not the same as just unique
    if do_name_replace
        unique_names_str = [string(k, base=2) for k in unique_names]
        unique_names_str = [repeat("0", longest_name - length(c)) * c for c in unique_names_str]
        if unique_names_str == clusters_names
            for (fixed_name, cname) in zip(unique_names_str, clusters_names)
                if fixed_name != cname
                    clust = pop!(clusters_merges, cname)

                    clusters_merges[fixed_name] = ClusterSplit(fixed_name, clust.range, clust.deph, clust.height)
                end
            end
        end
    end


    # is correction needed?
    if findmin(clusters_names_len)[1] != findmax(clusters_names_len)[1]
        target_length = findmax(clusters_names_len)[1]

        for k in clusters_names
            name_len = k |> length
            if name_len != target_length
                clust = pop!(clusters_merges, k)
                fixed_name = clust.name
                for n in 1:(target_length-name_len)
                    fixed_name *= "0"
                end
                clusters_merges[fixed_name] = ClusterSplit(fixed_name, clust.range, clust.deph, clust.height)
            end
        end
    end # if different

    clusters_merges = clusters_merges |> sort
    return clusters_merges
end


"""
Get cluster merges from cutree.

This function is not finished
"""
function get_cluster_merges(extended_cluster, cut_height)

    unique_clusters = merged_cluster[selected_cluster.order] |> unique
    inv_leaves_range_per_merge = extended_cluster.ranges[end:-1:1]


    first_used = findfirst(x -> x == "1", ["$(k)" for k in bitstring(findmax(unique_clusters)[1])])
    used_range = max(first_used - 1, 1)


    clusters_merges = Dict()

    cluster_name = unique_clusters[1]
    for cluster_name in unique_clusters

        clust_name = bitstring(cluster_name)[used_range:end]

        # Find range for each cluster
        cluster_range = inv_leaves_range_per_merge[1]

        ClusterSplit(clust_name, cluster_range, 0, cut_height)
    end


    return clusters_merges
end

"""
    translate_merges_to_treecut(clusters_merges, selected_cluster::Hclust)

Returns the the vector describing membership of cluster leaves after merging
(from top) n-branches of the cluster

Returned value has the same structure as if the `cuttree` was applied to
a `Hclust` instance.
"""
function translate_merges_to_cutree(clusters_merges, selected_cluster::Hclust)
    merging = zeros(Int, length(selected_cluster.order))
    for (k, c) in clusters_merges
        cluster_position = c.range |> collect
        cluster_id = get_cluster_split_name(c)
        for position in cluster_position
            merging[position] = cluster_id
        end
    end
    return merging
end


function chceck_yeo_representation(input_matrix::Matrix{T}, brain_regions_info) where {T}
    total_row, total_cols = size(input_matrix)
    yeo_representation = zeros(T, 8, total_cols)

    yeo_vector = [v[2] for (k, v) in brain_regions_info]
    yeo_networks = yeo_vector |> unique
    network_positions = [net => findall(x -> x == net, yeo_vector) for net in yeo_networks] |> Dict

    @debug (vcat([v for (k, v) in network_positions]...) |> sort == 1:total_row)

    for subject in 1:total_cols
        for (k, net) in yeo_networks |> enumerate
            val = 0
            if (@pipe net |> network_positions[_] |> input_matrix[_, subject] |> sum) > 0
                val = 1
            else
                val = 0
            end
            yeo_representation[k, subject] = val
        end
    end

    return yeo_representation
end
