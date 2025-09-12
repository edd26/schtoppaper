
"""
    plot_coloured_cluster_merge(cluster, height::Int=5)

Return a cluster plot with background coloured according to the `cluster` merge at level `height`.
Source: https://discourse.julialang.org/t/statsplots-dendrogram-how-to-have-clustered-coloring/50768/12
"""
function plot_coloured_cluster_merge(selected_cluster; height::Int=5)
    y_hat = cutree(selected_cluster, h=height)

    extremums = Vector(undef, length(unique(y_hat)))
    for (k, i) in enumerate(unique(y_hat))
        extremums[k] = extrema(
            [findfirst(x -> x .== j, selected_cluster.order) for j in findall(x -> x == i, y_hat)]
        ) .+ (-0.5, 0.5)
    end
    coloured_clust_plt = plot([Shape([m, M, M, m], [0, 0, height, height]) for (m, M) in extremums];
        alpha=0.5
    )
    plot!(selected_cluster,
        xticks=(1:length(y_hat),
            string.(selected_cluster.order)),
        xlabel="node",
        ylabel=string(selected_cluster.linkage) * " linkage distance"
    )

    return coloured_clust_plt
end


function plot_coloured_cluster_merge2(selected_cluster; min_clusters::Int=5)
    y_hat = cutree(selected_cluster, k=min_clusters)

    extremums = Vector(undef, length(unique(y_hat)))
    for (k, i) in enumerate(unique(y_hat))
        extremums[k] = extrema(
            [findfirst(x -> x .== j, selected_cluster.order) for j in findall(x -> x == i, y_hat)]
        ) .+ (-0.5, 0.5)
    end
    coloured_clust_plt = plot([Shape([m, M, M, m], [0, 0, 20, 20]) for (m, M) in extremums];
        alpha=0.5
    )
    plot!(selected_cluster,
        xticks=(1:length(y_hat),
            string.(hcl1.order)),
        xlabel="node",
        ylabel=string(hcl1.linkage) * " linkage distance"
    )

    return coloured_clust_plt
end

"""

Plots cluster (with colours indicating merges) 
"""
function get_cluster_with_colouring(selected_cluster, merged_cluster, selected_palette, cluster_height, min_clusters)
    ## ===-===-===-===-===-
    # Cluster plot
    extremums = [
        extrema([findfirst(x -> x .== j, selected_cluster.order) for j in findall(x -> x == i, merged_cluster)]
        ) .+ (-0.5, 0.5)
        for (k, i) in enumerate(unique(merged_cluster))]

    matrix_ranges = [ceil(Int, m):ceil(Int, M - 1) for (m, M) in extremums] |> sort
    total_ranges = max(length(matrix_ranges), 2)

    palette_args = (selected_palette, total_ranges)
    palette_kwargs = (rev=true,)
    max_height = max(selected_cluster.height...)
    ystep = min(0.2 * max_height, 20)

    # ===-===-===-===-===-
    # Cluster plot
    coloured_clust_plt = Plots.plot(
        [Shape([m, M, M, m], [0, 0, cluster_height, cluster_height]) for (m, M) in extremums];
        alpha=0.5,
        palette=palette(palette_args...; palette_kwargs...),
        legend=:outerright,
        legend_column=total_ranges,
        labels=""
    )
    Plots.plot!(selected_cluster,
        xticks=false,
        yticks=0:ystep:max_height,
        ylabel=string(selected_cluster.linkage) * " linkage distance",
        title="linkage= $(selected_cluster.linkage), max height=$(cluster_height) min clusters=$(min_clusters)",
        legend=:outerright,
        cbar=false,
        labels="",
    )
    return coloured_clust_plt, extremums, matrix_ranges
end

"""

Plots cluster (with colours indicating merges) with region matrix indicatig brin region
presence in the cycles.
"""
function get_cluster_with_regions_distribution(
    selected_cluster,
    regions_in_cycles_matrix,
    selected_linkage,
    unique_clusters;
    final_plt_proportions::Vector{Float64}=[0.4, 0.6],
    selected_palette=:matter
)
    coloured_clust_plt, extremums, matrix_ranges = get_cluster_with_colouring(selected_cluster, merged_cluster, selected_palette, cluster_height, min_clusters)
    total_ranges = length(matrix_ranges)
    palette_args = (selected_palette, total_ranges)
    palette_kwargs = (rev=true,)

    # ===-===-===-===-===-===-===-===-===-===-===-===-
    # Raster plot preparation
    regions_vs_cycles_for_large_scale = large_scale_network_matrix_preprocessing(regions_in_cycles_matrix, yeo7_networks, brain_regions_info)

    xs1 = range(1, stop=size(regions_vs_cycles_for_large_scale, 2))
    ys1 = ["$(new_regions_order[y])" for y in 1:total_regions]
    xtick_step = 20

    x_labels = vcat([["y$(id)" for r in rang] for (id, rang) in zip(unique_clusters, matrix_ranges)]...)

    # ===-===-===-===-===-
    # Raster plot plotting
    hmp_plt = Plots.heatmap(xs1,
        ys1,
        regions_vs_cycles_for_large_scale[new_regions_order, selected_cluster.order];
        alpha=0.1,
        c=:white,
        cbar=false
    )
    Plots.plot!(
        [Shape([m, M, M, m], [0, 0, 94, 94]) for (m, M) in extremums];
        alpha=0.5,
        palette=palette(selected_palette, total_ranges + 1; palette_kwargs...),
        legend=:outerright,
        legend_column=total_ranges
    )
    hmp_plt = Plots.heatmap!(xs1,
        ys1,
        regions_vs_cycles_for_large_scale[new_regions_order, selected_cluster.order];
        c=colours_palete,
        cbar=false,
        dpi=300,
        ylabel="Region index",
        widen=false,
        xticks=(1:xtick_step:length(merged_cluster),
            x_labels[1:xtick_step:end]
        ),
        xlabel="node",
        xrotation=45,
        labels="",
        legend=:outerright
    )

    # ===-===-===-===-===-
    # Cluster plot and raster plot
    fin_plt = Plots.plot(
        coloured_clust_plt,
        hmp_plt;
        layout=grid(2, 1, heights=final_plt_proportions),
        size=(1200, 700),
        dpi=300,
        widen=false,
        thickness_scaling=1.4,
        link=:both
    )
    return fin_plt
end
