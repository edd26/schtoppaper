import Printf: @sprintf
using CairoMakie
using Makie.GeometryBasics
using PersistenceLandscapes
using DataStructures: OrderedDict
using Pipe

"ClusterMerging.jl" |> srcdir |> include

# ===-===-===-
# Makie functions for plotting
function add_cluster_background_v1!(ax_regions, cluster_ranges, color_vals; map_colours=:grays)
    map1 = NaN

    total_unique_clusters = cluster_ranges |> length
    for (color_id, (clust_key, clust_range),) in cluster_ranges |> enumerate

        left, right = clust_range
        map1 = makie_backend.vspan!(
            ax_regions,
            left,
            right,
            color=(color_vals[color_id], 0.4),
            transparency=true,
            # following two are required for legend creation
            colormap=Reverse(:grays),
            colorrange=(1, total_unique_clusters),
        )
    end
    return map1
end


function add_cluster_background!(ax_regions::Axis, merged_clusters::OrderedDict,
    cluster_info::OrderedDict; alpha=0.4, map_colours=Reverse(:grays), use_bcg=true, makie_backend=CairoMakie)
    map1 = NaN
    clust_names = [n for (n, v) in merged_clusters]

    max_cluster_index = clust_names |> length
    for cname in clust_names
        r = merged_clusters[cname].range
        left, right = r[1], r[end]
        if left == right
            left -= 0.5
            right += 0.5
        end
        if use_bcg
            c = cluster_info[cname].background_colour
        else
            c = cluster_info[cname].signifficance_colour
        end

        map1 = makie_backend.vspan!(
            ax_regions,
            left,
            right,
            color=(c, alpha),
            transparency=true,
        )
    end
    return map1
end

function add_dendrogram_background!(ax_regions::Axis, merged_clusters::OrderedDict, cluster_info::OrderedDict; alpha=0.4, use_signifficance_colour=false)
    map1 = NaN

    clust_names = [n for (n, v) in merged_clusters]
    for cname in clust_names
        r = merged_clusters[cname].range
        clust_height = merged_clusters[cname].height
        left, right = r[1], r[end]
        if use_signifficance_colour
            bcg_colour = cluster_info[cname].signifficance_colour
        else
            bcg_colour = cluster_info[cname].background_colour
        end

        p = Polygon(Point2f[(left, 0), (right, 0), (right, clust_height), (left, clust_height)],)
        map1 = poly!(
            ax_regions,
            p,
            color=(bcg_colour, alpha),
            transparency=true,
        )
    end
    return map1
end


function add_dendrogram_background!(ax_regions, merged_cluster, clusters_index_ranges, color_vals, height)
    map1 = NaN

    total_unique_clusters = (merged_cluster) |> unique |> length
    for (color_id, clust_range,) in clusters_index_ranges
        left, right = clust_range .- 0.5
        p = Polygon(Point2f[(left, 0), (right, 0), (right, height), (left, height)],)
        map1 = poly!(p,
            color=(color_vals[color_id], 0.4),
            transparency=true,
        )
    end
    return map1
end



function add_dendrogram_background!(ax_regions, hclust_m; col_palette=[])
    map1 = NaN
    max_cluster_index = get_max_cluster_index(hclust_m)
    if col_palette |> isempty
        col_palette = palette(:grays, max_cluster_index; rev=true)
    end

    for (k, (name, clust)) in hclust_m.clusters |> enumerate

        left, right = (clust.range[1], clust.range[end]) .- 0.5
        left -= 0.25
        right += 0.25

        p = Polygon(Point2f[(left, 0), (right, 0), (right, clust.height), (left, clust.height)],)
        map1 = poly!(ax_regions, p, color=(col_palette[k], 0.4), transparency=true,)
    end
    return map1
end

function get_max_cluster_index(hclust_m)
    max_cluster_index = findmax([parse(Int, k, base=2) for (k, v) in hclust_m.clusters])[1] + 1
    return max_cluster_index
end


# https://github.com/MakieOrg/Makie.jl/issues/398
function treepositions(hc::Union{Hclust,HclustMerged}; useheight=true, orientation=:vertical, kwargs...)
    order = StatsBase.indexmap(hc.order)
    nodepos = Dict(-i => (float(order[i]), 0.0) for i in hc.order)
    xs = []
    ys = []
    for i in 1:size(hc.merges, 1)
        x1, y1 = nodepos[hc.merges[i, 1]]
        x2, y2 = nodepos[hc.merges[i, 2]]
        xpos = (x1 + x2) / 2
        ypos = useheight ? hc.heights[i] : (max(y1, y2) + 1)
        nodepos[i] = (xpos, ypos)
        push!(xs, [x1, x1, x2, x2])
        push!(ys, [y1, ypos, ypos, y2])
    end
    if orientation == :horizontal
        return ys, xs
    else
        return xs, ys
    end
end

function dendrogram!(axis, h; color=:blue, kwargs...)
    CairoMakie.Axis(Figure()[1, 1])
    for (x, y) in zip(treepositions(h;)...)
        lines!(axis, x, y; color, kwargs...)
    end
end

function look_up_colour(x, extremas)
    for (id, (start_x, stop_x)) in extremas |> enumerate
        if x > start_x && x < stop_x
            return id
        else
            ErrorException("Falied to determain region")
        end
    end

end

function dendrogram_coloured!(axis, h; color=:blue, kwargs...)
    for (id, (x, y)) in zip(treepositions(h;)...) |> enumerate
        color = :blue
        lines!(axis, x, y; color=color, kwargs...)
    end
end

function dendrogram_coloured!(axis, h, max_height, cluster_ranges, col_palette; color=:blue, kwargs...)
    keys_ = [k for k in cluster_ranges |> keys] |> sort
    extremas = [cluster_ranges[k] for k in keys_]
    for (id, (x, y)) in zip(treepositions(h;)...) |> enumerate
        if y[2] < max_height
            col_id = look_up_colour(x[1], extremas)
            color = col_palette[col_id]
            lines!(axis, x .- 0.5, y; color, kwargs...)
        else
            color = :blue
            lines!(axis, x, y; color, kwargs...)
        end
    end
end

function dendrogram_coloured!(axis, hclust_m::HclustMerged; col_palette=[], skip_naming=false, kwargs...)
    total_clusters = get_max_cluster_index(hclust_m)

    if col_palette |> isempty
        col_palette1 = palette(:brocO, total_clusters;)
        col_palette = []
        for k in 1:(total_clusters÷2)
            push!(col_palette, col_palette1[k])
            push!(col_palette, col_palette1[end-k+1])
        end
        col_palette = [k for k in col_palette]
        col_palette = palette(:grays, total_clusters; rev=true)
    end

    tree_positions = treepositions(hclust_m;)
    for (x, y) in zip(tree_positions...)

        leaves_range = ceil(Int, x[1]):ceil(Int, x[3])
        if skip_naming
            cluster_name = ""
        else
            cluster_name = find_superset_from_clusters(leaves_range, hclust_m.clusters)
        end

        if cluster_name == ""
            color = :red
            lines!(axis, x, y; color, kwargs...)
        else
            if all(y .== 0)
                local drawing_x = x .+ 1
            else
                local drawing_x = x
            end
            related_clust = hclust_m.clusters[cluster_name]
            if y[2] <= related_clust.height
                col_id = parse(Int, related_clust.name, base=2) + 1
                color = col_palette[col_id]
                lines!(axis, drawing_x .- 0.5, y; color, kwargs...)
            else
                color = :blue
                lines!(axis, drawing_x, y; color, kwargs...)
            end
        end
    end
end

function plot_persistence_landscape!(
    plt_axis,
    pl1::PersistenceLandscape;
    max_layers=size(pl1.land, 1),
    max_colour_range=size(pl1.land, 1),
    alpha::Float64=0.0,
    colors=nothing,
    plot_kwargs...
)
    if max_colour_range < max_layers
        @warn "Selected colour range is less than total layers! Changing to max layers instead"
        max_colour_range = max_layers
    end

    colors =
        if isnothing(colors)
            [RGBf(c) for c in cgrad(:cmyk, max(2, max_colour_range), categorical=true, rev=true, alpha=alpha)]
        else
            colors
        end

    try
        colors = cgrad(plot_kwargs[:palette], max_colour_range, categorical=true, rev=true, alpha=alpha)
    catch
        @debug "Catched no palette"
    end

    for k = 1:max_layers
        peaks_position, peaks = PersistenceLandscapes.get_peaks_and_positions(pl1.land[k])

        CairoMakie.lines!(
            plt_axis,
            peaks_position,
            peaks;
            color=colors[k],
            plot_kwargs...
        )
    end
end

function mplot_vertical_barcodes!(ax_barcodes, bd_data, cycles_clustering_order; alpha_lvl=0.1, color_bank=[:blue, :orange, :red])
    # Prepare barcodes
    used_keys = bd_data |> keys |> collect

    bd_matrix_dim1 = bd_data[used_keys[1]]
    if length(used_keys) > 1
        local_key = used_keys[2]
        for local_key in used_keys[2:end]
            bd_matrix_dim1 = hcat(bd_matrix_dim1, bd_data[local_key])
        end
    end

    # Do the plot
    for (index, cycle_pointer) in enumerate(cycles_clustering_order)

        bd_vector = bd_matrix_dim1[cycle_pointer, :]
        for (k, persistence_key) in used_keys |> enumerate

            for bd_pair in bd_vector[k]
                if occursin("hc", persistence_key)
                    bar_color = color_bank[1]
                elseif occursin("sch", persistence_key)
                    bar_color = color_bank[2]
                else
                    @warn "Unknown persistence key"
                    bar_color = color_bank[k]
                end

                vlines!(
                    ax_barcodes,
                    index,
                    ymin=bd_pair[1],
                    ymax=bd_pair[2],
                    color=(bar_color, alpha_lvl),
                    label="",
                    alpha=alpha_lvl,
                    transparency=true
                ) # same low and high error

            end
        end
    end

end


function get_landscapes_and_histograms_plots!(
    gc,
    distances_distribution,
    landscapes_distance,
    total_cols,
    total_rows,
    nth_percentile,
    total_bins,
    cluster_plot_info,
    subject_related_barcodes,
    average_landscape_per_cluster,
    data_keys,
    hclust_m;
    do_two_tailed=true,
    plt_height=100,
    plt_width=400,
    sort_by_position=true
)
    hist, vline1, vline2, vline3 = 0, 0, 0, 0

    clust_colouring = ["clust$(c |> get_cluster_split_name)" => c.signifficance_colour
                       for (cname, c) in cluster_plot_info
    ] |> OrderedDict

    if do_two_tailed
        left_margin = @sprintf "%0.2f" 100 * (1 - nth_percentile) / 2
        right_margin = @sprintf "%0.2f" 100 * (nth_percentile + (1 - nth_percentile) / 2)
    else
        right_margin = @sprintf "%0.2f" nth_percentile * 100
        left_margin = NaN
    end

    plt_col_index = 1
    plt_row_index = 1
    all_keys = [k for k in keys(distances_distribution)]
    if sort_by_position
        ordered_permutation = [mean(c.range) for (k, c) in hclust_m.clusters] |> sortperm
    else # sort by number
        keys_numbers = [parse(Int, split(k, "clust")[2]) for k in all_keys]
        ordered_permutation = sortperm(keys_numbers)
    end
    for (cl_id, cluster_key) in all_keys[ordered_permutation] |> enumerate
        backgroundcolor_val = clust_colouring[cluster_key]

        # Get the statistical information
        left_val, right_val = get_nth_percentile(
            distances_distribution[cluster_key],
            NTH_PERCENTILE;
            two_tailed=DO_TWO_TAILED
        )

        # set up row and col for plt
        gbcurrent_plt = gc[plt_row_index, plt_col_index] = GridLayout() # cluster and raster
        plt_col_index += 1
        if plt_col_index > total_cols
            plt_col_index = 1
            plt_row_index += 1
        end
        if plt_row_index > total_rows
            @warn "Unable to fit plots in the space"
        end

        # ===-===-
        # get landscape plots
        if ismissing(plt_width) || ismissing(plt_height)
            ax_land = [
                Axis(gbcurrent_plt[1, 1], backgroundcolor=backgroundcolor_val),
                Axis(gbcurrent_plt[2, 1], backgroundcolor=backgroundcolor_val)
            ]
            ax_hist_plt = Axis(gbcurrent_plt[3, 1])
        else
            ax_land = [
                Axis(gbcurrent_plt[1, 1], backgroundcolor=backgroundcolor_val, width=plt_width, height=plt_height),
                Axis(gbcurrent_plt[2, 1], backgroundcolor=backgroundcolor_val, width=plt_width, height=plt_height)
            ]
            ax_hist_plt = Axis(gbcurrent_plt[3, 1], width=plt_width, height=plt_height,)
        end

        for (id, key_related) in data_keys |> enumerate
            key_related_barcodes = [bar for (k, bar) in subject_related_barcodes[cluster_key][key_related] if length(bar) > 0]
            total_key_related_members = key_related_barcodes |> length
            total_key_related_bars = vcat(key_related_barcodes...) |> length
            key_related_title = "$(key_related), $(cluster_key), ppl=$(total_key_related_members), bars=$(total_key_related_bars)"

            plot_persistence_landscape!(
                ax_land[id],
                average_landscape_per_cluster[cluster_key][key_related],
                max_colour_range=9,
                linewidth=3,
            )
            ax_land[id].title = key_related_title
            makie_backend.xlims!(ax_land[id], 0, x_max)
            makie_backend.ylims!(ax_land[id], 0, y_max)


            id == 1 && hidexdecorations!(ax_land[id])#, grid=false)
        end

        # ===-===-
        # get histogram plot
        hist = makie_backend.hist!(
            ax_hist_plt,
            distances_distribution[cluster_key],
            bins=total_bins,
            normalization=:pdf,
            label="distance distribution",
        )
        vline1 = makie_backend.vlines!(
            ax_hist_plt,
            [landscapes_distance[cluster_key]],
            color=(:red),
            label="sch-hc",
            linewidth=6,
        )
        if left_margin |> isnan
            vline2 = Nothing
        else
            vline2 = makie_backend.vlines!(
                ax_hist_plt,
                [left_val],
                color=:orange,
                label="$(left_margin)% of data",
                linewidth=4,
            )
        end
        vline3 = makie_backend.vlines!(
            ax_hist_plt,
            [right_val],
            color=:orange,
            label="$(right_margin)% of data",
            linewidth=4,
        )
        ax_hist_plt.xlabel = "landscape distance"
        ax_hist_plt.ylabel = "pdf"

        Makie.xlims!(
            ax_hist_plt,
            low=0,
        )
        Makie.ylims!(
            ax_hist_plt,
            low=0,
        )

    end

    return hist, vline1, vline2, vline3
end



function visualise_dendrogram(data_clust, merged_cluster;
    final_plt_width=900,
    final_plt_height=300)
    f = Figure(
        resolution=(final_plt_width, final_plt_height)
    )

    # Setting up GridLayouts
    fgl = f[1, 1] = GridLayout() # whole space
    ga = fgl[1, 1] = GridLayout() # main plots

    # Panel A - cluster and raster
    # Find colour ranges
    cluster_ranges = [
        extrema([findfirst(x -> x .== j, data_clust.order) for j in findall(x -> x == i, merged_cluster)]
        ) .+ (-0.5, 0.5)
        for (k, i) in enumerate(unique(merged_cluster))]


    unique_clusters = merged_cluster[data_clust.order] |> unique
    total_clusters = length(unique_clusters)
    col_palette = palette(:grays, total_clusters; rev=true)
    color_vals = [RGBf(k) for k in col_palette]

    # ===-
    # Plot dendrogram
    ax_dendrogram = Axis(ga[1, 1])

    dendrogram_coloured!(
        ax_dendrogram,
        data_clust,
        100,
        cluster_ranges,
        col_palette;
        color=:blue
    )

    add_dendrogram_background!(ax_dendrogram, merged_cluster, cluster_ranges, color_vals, cluster_height)

    # Set up parameters for dendrogram plot 
    ax_dendrogram.xticks = 0:0
    ax_dendrogram.xticklabelrotation = pi / 4
    max_height = max(data_clust.height...)
    makie_backend.ylims!(ax_dendrogram, 0, 1.1max_height)
    makie_backend.xlims!(ax_dendrogram, 0, length(merged_cluster))
    ax_dendrogram.yticks = 0:50:1.1max_height

    hidexdecorations!(ax_dendrogram, grid=false)

    return f
end


"""
    get_density_estimation_plot!(axis, clust_info, estimation_data)

Plots a series of density estimation plots, one for each unique element in 'clust_info'.
Results are similar to violin plots.
"""
function get_density_estimation_plot!(axis, clust_info, estimation_data)
    for clust in clust_info |> unique
        plt_data_position = clust_info .== clust
        plt_data = estimation_data[plt_data_position]
        makie_backend.density!(
            axis,
            plt_data,
            direction=:y,
            offset=clust,
            npoints=plt_data |> length,
            color=:y,
            strokewidth=1,
            strokecolor=:black
        )
    end
end



"""
get_background_colouring(landscapes_distance, distances_distribution_sorted)

Produces backgrounds for each of the cluster. This is done by getting the FDR test on each of the landscapes
    distribution from `distances_distribution`.

Possible backgrounds are:
- red for clusters for which the FDR test was passed;
- gray fro clusters not passing FDR but having low p-value;
- white otherwise.

Returns dictionary with keys being clusters numbers and values corresponding backgrounds.
"""
function get_background_colouring(landscapes_distance, distances_distribution_sorted; intense_signifficance=false)
    dendrogram_sorted_keys = [k for (k, v) in landscapes_distance]

    critical_values = get_BH_critical_values_form_landscapes(landscapes_distance, distances_distribution_sorted)
    subselected_cluster_colours = OrderedDict()

    for ck in dendrogram_sorted_keys
        ck_index = parse(Int, split(ck, "clust")[2])

        cluster_position = critical_values.ClusterKeys .== ck
        interesting_row = critical_values[cluster_position, :]

        signifficant_check = interesting_row[!, "small p value"][1]
        BH_procedure_passed = interesting_row[!, :AcceptedValue][1]

        if BH_procedure_passed
            if intense_signifficance
                backgroundcolor_val = RGB(([185, 65, 72] ./ 255)...)
            else
                backgroundcolor_val = RGB(1.0, 0.8, 0.8)
            end

        elseif signifficant_check
            backgroundcolor_val = RGB(0.9, 0.9, 0.9)
        else
            backgroundcolor_val = RGB(1.00, 1.00, 1.00)
        end

        subselected_cluster_colours[ck_index] = backgroundcolor_val
    end
    return subselected_cluster_colours
end

struct ClusterPlotInfo
    name
    background_colour
    signifficance_colour
end

function get_cluster_split_name(clust_split::ClusterPlotInfo)
    return parse(Int, clust_split.name, base=2)
end
