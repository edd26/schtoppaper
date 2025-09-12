using Plots
using Plots.PlotMeasures
gr()

function large_scale_network_matrix_preprocessing(regions_in_cycles_matrix, yeo7_networks, brain_regions_info)
    regions_vs_cycles_for_large_scale = copy(regions_in_cycles_matrix)
    all_missing = findall(x -> ismissing(x), regions_vs_cycles_for_large_scale)
    regions_vs_cycles_for_large_scale[all_missing] .= 0
    non_0 = findall(x -> x > 0, regions_vs_cycles_for_large_scale)
    regions_vs_cycles_for_large_scale[non_0] .= 1

    large_scale_to_number = Dict()
    for (region_index, n_region) = enumerate(yeo7_networks)
        large_scale_to_number[n_region] = region_index
    end

    for (region_index, n_region) = brain_regions_info
        regions_vs_cycles_for_large_scale[region_index, :] .= large_scale_to_number[n_region[2]]
    end

    regions_vs_cycles_for_large_scale[all_missing] .= missing
    return regions_vs_cycles_for_large_scale
end

function large_networks_plotting(cycles_clustering,
    yeo7_networks,
    brain_regions_info,
    regions_vs_cycles_for_large_scale,
    new_regions_order,
    ordered_bd_type_split,
    regions_coloring;
    ppl_matrix=[],
    plt_size=(3508, 2480),
    dim=1,
    title_appendix="",
    total_regions=94
)
    colours_palete = [brain_regions_info[findfirst(x -> x[2] == region, brain_regions_info)][4] for region in yeo7_networks]
    # ===-===-===-===-
    # Cluster
    cluster_hmap = plot(cycles_clustering,
        ylabel="a. u.",
        title="Dimension = $(dim)" * title_appendix,
        label="",
        xticks=false)

    # ===-===-===-===-===-===-===-===-===-
    # views on the brain
    brain_plot2 = plot_brain_projection(yeo7_networks, brain_regions_info, regions_coloring)

    # ===-===-===-===-===-===-===-===-===-===-
    # Raster plot
    xs1 = range(1, stop=size(regions_vs_cycles_for_large_scale, 2))
    ys1 = ["$(new_regions_order[y])" for y in 1:total_regions]
    hmp_plt = heatmap(xs1,
        ys1,
        regions_vs_cycles_for_large_scale[new_regions_order, cycles_clustering.order];
        c=colours_palete,
        cbar=false,
        dpi=300,
        xticks=0:50:size(regions_vs_cycles_for_large_scale, 2),
        ylabel="region index"
    )

    # ===-===-===-===-===-===-===-===-===-===-
    # Regions legend
    xs = ones(Float64, 1)
    ys = ["$(value)__$(brain_regions_info[value][1])" for value in new_regions_order]
    regions_legend = heatmap(xs,
        ys,
        regions_coloring[new_regions_order, :],
        ticks=:all,
        xticks=:false, # 0:1:1.5,
        colorbar=:false,
        c=colours_palete,
        tickfontsize=6,
        guide_position=:right,
        widen=false,
    )

    # ===-===-===-===-===-===-===-===-===-===-
    # Empty plot
    empty_plot = plot(ticks=nothing, border=:none)

    # ===-===-===-===-===-===-===-===-===-===-
    # vertical barcodes
    bd_canvas = plot_vertical_barcodes(ordered_bd_type_split, cycles_clustering.order)

    # # ===-===-===-===-===-===-===-===-===-===-
    # final plot
    final_kwargs = (
        size=plt_size,
        dpi=300,
        widen=false,
        thickness_scaling=1.2,
        left_margin=10PlotMeasures.mm,
        bottom_margin=10PlotMeasures.mm,
    )
    # #
    if !isempty(ppl_matrix)
        plt4, plt4_complementary = plot_ppl_matrix(ppl_matrix, cycles_clustering)

        return map_brain_plt = plot(cluster_hmap,
            brain_plot2,
            hmp_plt,
            regions_legend,
            plt4,
            plt4_complementary,
            bd_canvas,
            empty_plot;
            layout=grid(4, 2, heights=[0.3, 0.4, 0.2, 0.1], widths=[0.8, 0.2]),
            final_kwargs...
        )
    else
        # ===-===-===-===-===-===-===-===-===-===-
        return map_brain_plt = plot(cluster_hmap,
            brain_plot2,
            hmp_plt,
            regions_legend,
            bd_canvas,
            empty_plot;
            layout=grid(3, 2, heights=[0.3, 0.5, 0.2], widths=[0.8, 0.2]),
            final_kwargs...
        )
    end
end


function large_net_merging_comparison_plotting(cycles_clustering,
    yeo7_networks,
    brain_regions_info,
    regions_vs_cycles_for_large_scale,
    regions_vs_cycles_for_large_scale_original,
    new_regions_order,
    ordered_bd_type_split,
    regions_coloring;
    ppl_matrix=[],
    plt_size=(3508, 2480),
    dim=1,
    title_appendix=""
)

    colours_palete = [brain_regions_info[findfirst(x -> x[2] == region, brain_regions_info)][4] for region in yeo7_networks]
    # ===-===-===-===-
    # Cluster
    cluster_hmap = plot(cycles_clustering,
        ylabel="a. u.",
        title="Dimension = $(dim)" * title_appendix,
        label="",
        xticks=false)

    # ===-===-===-===-===-===-===-===-===-
    # views on the brain
    brain_plot2 = plot_brain_projection(yeo7_networks, brain_regions_info, regions_coloring)

    # ===-===-===-===-===-===-===-===-===-===-
    # Clustered raster plot
    xs1 = range(1, stop=size(regions_vs_cycles_for_large_scale, 2))
    ys1 = ["$(new_regions_order[y])" for y in 1:total_regions]
    hmp_plt = heatmap(xs1,
        ys1,
        regions_vs_cycles_for_large_scale[new_regions_order, cycles_clustering.order];
        c=colours_palete,
        cbar=false,
        dpi=300,
        xticks=0:50:size(regions_vs_cycles_for_large_scale, 2),
        ylabel="region index"
    )

    # ===-===-===-===-===-===-===-===-===-===-
    # Original raster plot
    xs2 = range(1, stop=size(regions_vs_cycles_for_large_scale_original, 2))
    ys2 = ["$(new_regions_order[y])" for y in 1:total_regions]
    hmp_plt_original = heatmap(xs2,
        ys2,
        regions_vs_cycles_for_large_scale_original[new_regions_order, cycles_clustering.order];
        c=colours_palete,
        cbar=false,
        dpi=300,
        xticks=0:50:size(regions_vs_cycles_for_large_scale_original, 2),
        ylabel="region index",
        titla="Original regions vs cycles"
    )

    # ===-===-===-===-===-===-===-===-===-===-
    # Regions legend
    xs = ones(Float64, 1)
    ys = ["$(k)__$(v[1])" for (k, v) in brain_regions_info]
    regions_legend = heatmap(xs,
        ys,
        regions_coloring[new_regions_order, :],
        ticks=:all,
        xticks=:false, # 0:1:1.5,
        colorbar=:false,
        c=colours_palete,
        tickfontsize=6,
        guide_position=:right,
        widen=false,
    )

    # ===-===-===-===-===-===-===-===-===-===-
    # Empty plot
    empty_plot = plot(ticks=nothing, border=:none)

    # ===-===-===-===-===-===-===-===-===-===-
    # vertical barcodes
    bd_canvas = plot_vertical_barcodes(ordered_bd_type_split, cycles_clustering.order)

    # # ===-===-===-===-===-===-===-===-===-===-
    # final plot
    final_kwargs = (
        size=plt_size,
        dpi=300,
        widen=false,
        thickness_scaling=1.2,
        left_margin=10PlotMeasures.mm,
        bottom_margin=10PlotMeasures.mm,
    )
    if !isempty(ppl_matrix)
        plt4, plt4_complementary = plot_ppl_matrix(ppl_matrix, cycles_clustering)

        return map_brain_plt = plot(cluster_hmap,
            brain_plot2,
            hmp_plt,
            regions_legend,
            hmp_plt_original,
            empty_plot,
            plt4,
            plt4_complementary,
            bd_canvas,
            empty_plot;
            layout=grid(5, 2, heights=[0.2, 0.25, 0.25, 0.2, 0.1], widths=[0.8, 0.2]),
            final_kwargs...
        )
    else
        # ===-===-===-===-===-===-===-===-===-===-
        return map_brain_plt = plot(cluster_hmap,
            brain_plot2,
            hmp_plt,
            regions_legend,
            hmp_plt_original,
            empty_plot,
            bd_canvas,
            empty_plot;
            layout=grid(4, 2, heights=[0.2, 0.2, 0.2, 0.2], widths=[0.8, 0.2]),
            final_kwargs...
        )
    end
end

function plot_vertical_barcodes(ordered_persistence_data, cycles_clustering_order::Vector; kwargs...)
    used_keys = collect(keys(ordered_persistence_data))
    bd_matrix_dim1 = ordered_persistence_data[used_keys[1]]
    if length(used_keys) > 1
        for key in used_keys[2:end]
            bd_matrix_dim1 = hcat(bd_matrix_dim1, ordered_persistence_data[key])
        end
    end
    total_keys = length(used_keys)

    bd_canvas = Plots.plot(;
        xticks=0:50:size(bd_matrix_dim1, 1),
        xlims=(0, size(bd_matrix_dim1, 1)),
        widen=false,
        title="Vertical barcodes",
        xlabel="Cycle index",
        ylabel="Normalized lifetime",
        kwargs...
    )

    for (index, cycle_pointer) in enumerate(cycles_clustering_order)
        x_vals = ones(Int, 2) .* index

        bd_vector = bd_matrix_dim1[cycle_pointer, :]
        for (k, persistence_key) in used_keys |> enumerate
            for bd_pair in bd_vector[k]
                if persistence_key == "hc"
                    bar_color = :blue
                elseif persistence_key == "sch"
                    bar_color = :orange
                else
                    @warn "Unknown persistence key"
                    bar_color = :red
                end
                Plots.plot!(x_vals, bd_pair,
                    label="",
                    markeralpha=0.1,
                    alpha=0.1,
                    color=bar_color)
            end
        end
    end
    return bd_canvas
end

function plot_ppl_matrix(ppl_matrix, cycles_clustering)
    empty_plot = plot(ticks=nothing, border=:none)

    max_xtick = size(ppl_matrix, 2)
    max_ytick = size(ppl_matrix, 1)

    bar_data = copy(ppl_matrix)
    all_missing = findall(x -> ismissing(x), bar_data)
    if !isempty(all_missing)
        bar_data[all_missing] .= 0
    end
    all_2 = findall(x -> x == 2, bar_data)
    if !isempty(all_2)
        bar_data[all_2] .= 1
    end
    section = [1:max_ytick÷2, max_ytick÷2+1:max_ytick]

    part1 = [sum(bar_data[section[1], k]) for k in 1:max_xtick]
    part2 = [sum(bar_data[section[2], k]) for k in 1:max_xtick]
    unique_nloops = hcat(part1, part2)

    plt4 = plot(heatmap(ppl_matrix[:, cycles_clustering.order];
            c=palette([:blue, :orange], 2),
            yticks=(0:10:max_ytick),
            ylabel="Subject index",
            xlims=(1, max_xtick),
            legend=false,
            xticks=false
        ),
        groupedbar(unique_nloops[cycles_clustering.order, :];
            bar_position=:stack,
            dpi=300,
            xticks=false,
            xlims=(1, max_xtick),
            bar_width=1,
            yflip=true,
            lw=0,
            legend=false,
            top_margin=-12PlotMeasures.mm,
            border_style=:origin,
            ylabel="#"
        ),
        layout=grid(2, 1, heights=[0.7, 0.3], link=:x,),
    )
    plt4_complementary = plot(
        bar(sum(bar_data, dims=2);
            yticks=(0:10:max_ytick),
            ylims=(0, max_ytick + 1),
            orientation=:horizontal,
            lw=0,
            ytick=false,
            yaxis=false,
            left_margin=-20PlotMeasures.mm
        ),
        empty_plot;
        layout=grid(2, 1, heights=[0.97, 0.03]),
        c=:blue,
        orientation=:horizontal,
        lw=0,
        ytick=false,
        yaxis=false,
        legend=false,
        left_margin=-20PlotMeasures.mm)

    return plt4, plt4_complementary
end

function plot_ppl_matrix_COBRE(ppl_matrix, cycles_clustering_order::Vector)
    empty_plot = plot(ticks=nothing, border=:none)

    max_xtick = size(ppl_matrix, 2)
    max_ytick = size(ppl_matrix, 1)

    bar_data = copy(ppl_matrix)
    all_missing = findall(x -> ismissing(x), bar_data)
    bar_data[all_missing] .= 0
    all_2 = findall(x -> x == 2, bar_data)
    bar_data[all_2] .= 1
    section = [1:max_ytick]

    unique_nloops = [sum(bar_data[section[1], k]) for k in 1:max_xtick]


    plt4 = plot(heatmap(ppl_matrix[:, cycles_clustering_order];
            c=palette([:blue, :orange], 2),
            yticks=(0:10:max_ytick),
            ylabel="Time index",
            xlims=(1, max_xtick),
            legend=false,
            xticks=false
        ),
        groupedbar(unique_nloops[cycles_clustering_order, :];
            bar_position=:stack,
            dpi=300,
            xticks=false,
            xlims=(1, max_xtick),
            bar_width=1,
            yflip=true,
            lw=0,
            legend=false,
            top_margin=-12PlotMeasures.mm,
            border_style=:origin,
            ylabel="#"
        ),
        layout=grid(2, 1, heights=[0.7, 0.3], link=:x,),
    )
    plt4_complementary = plot(
        bar(sum(bar_data, dims=2);
            yticks=(0:10:max_ytick),
            ylims=(0, max_ytick + 1),
            orientation=:horizontal,
            lw=0,
            ytick=false,
            yaxis=false,
            left_margin=-20PlotMeasures.mm
        ),
        empty_plot;
        layout=grid(2, 1, heights=[0.97, 0.03]),
        c=:blue,
        orientation=:horizontal,
        lw=0,
        ytick=false,
        yaxis=false,
        legend=false,
        left_margin=-20PlotMeasures.mm)

    return plt4, plt4_complementary
end

function plot_brain_projection(yeo7_networks, brain_regions_info, regions_coloring)
    x_label = "<- left  right ->"
    y_label = "<- ant.  post. ->"
    z_label = "<- sup.  inf. ->"

    brain_veiw_kwargs = (markeralpha=0.9,
        markerstrokewidth=0,
    )

    # ===-===-
    transverse_view = plot(title="Brain view, transverse plane",
        legend=false,
        xlabel=y_label,
        ylabel=x_label,
        aspect_ratio=1,
    )

    for (region_index, n_region) = enumerate(yeo7_networks)
        regions_related_coords = [value[3] for (key, value) in brain_regions_info if value[2] == n_region]

        large_region_relatives = [key for (key, value) in brain_regions_info if value[2] == n_region]
        regions_coloring[large_region_relatives, 1] .= region_index

        related_points = length(regions_related_coords)
        related_x = [regions_related_coords[k][1] for k in 1:related_points]
        related_y = [regions_related_coords[k][2] for k in 1:related_points]

        scatter!(related_y, related_x;
            c=colours_palete[region_index],
            label="",
            brain_veiw_kwargs...)
    end
    # ===-===-
    sagittal_view = plot(title="Brain view, sagittal plane",
        legend=false,
        xlabel=y_label,
        ylabel=z_label,
        aspect_ratio=1,
    )
    for (region_index, n_region) = enumerate(yeo7_networks)
        regions_related_coords = [value[3] for (key, value) in brain_regions_info if value[2] == n_region]

        large_region_relatives = [key for (key, value) in brain_regions_info if value[2] == n_region]
        regions_coloring[large_region_relatives, 1] .= region_index

        related_points = length(regions_related_coords)
        related_x = [regions_related_coords[k][1] for k in 1:related_points]
        related_y = [regions_related_coords[k][2] for k in 1:related_points]
        related_z = [regions_related_coords[k][3] for k in 1:related_points]

        scatter!(related_y, related_z;
            c=colours_palete[region_index],
            label="",
            brain_veiw_kwargs...)
    end
    # ===-===-
    coronal_view = plot(title="Brain view, coronal plane",
        legend=false,
        xlabel=x_label,
        ylabel=z_label,
        aspect_ratio=1,
    )
    for (region_index, n_region) = enumerate(yeo7_networks)
        regions_related_coords = [value[3] for (key, value) in brain_regions_info if value[2] == n_region]

        large_region_relatives = [key for (key, value) in brain_regions_info if value[2] == n_region]
        regions_coloring[large_region_relatives, 1] .= region_index

        related_points = length(regions_related_coords)
        related_x = [regions_related_coords[k][1] for k in 1:related_points]
        related_y = [regions_related_coords[k][2] for k in 1:related_points]
        related_z = [regions_related_coords[k][3] for k in 1:related_points]

        scatter!(related_x, related_z;
            c=colours_palete[region_index],
            label="",
            brain_veiw_kwargs...)
    end
    # ===-===-
    empty_legend_plane = plot(ticks=nothing, border=:none, legend=:bottomright, xlims=(0, 1), ylims=(0, 1))
    for (region_index, n_region) = enumerate(yeo7_networks)
        scatter!([0.0, 1.0], [0.0, 1.0]; c=colours_palete[region_index], label=n_region, brain_veiw_kwargs...)
    end
    # ===-===-
    brain_plot2 = plot(
        sagittal_view,
        coronal_view,
        transverse_view,
        empty_legend_plane,
        layout=grid(2, 2),
        left_margin=20PlotMeasures.mm,
        bottom_margin=10PlotMeasures.mm
    )
end
# Vertical bars legend
# orange_birth_r1 = 0.1
# orange_death_r1 = 0.2
# orange_birth_r2 = 0.2
# orange_death_r2 = 0.3
#
# blue_birth_r1 = 0.1
# blue_death_r1 = 0.2
# blue_birth_r2 = 0.2
# blue_death_r2 = 0.3
#
# blue_kwargs = (label="",
#                 alpha=0.2,
#                 lw=8,
#                 color=:blue)
#
#
# orange_kwargs = (label="",
#                 alpha=0.2,
#                 lw=8,
#                 color=:orange)
#
# barcodes_legend_canvas = plot(
#                 xticks=0:1:11,
#                 xlims=(0,11 ),
#                 widen=false,
#                 title="Vertical barcodes",
#                 # xlabel="Cycle index",
#                 ylabel="Normalized lifespan",
#                 )
# total_cols = 10
# half_total_cols = total_cols÷2
# for column = 1:total_cols
#     x_vals = ones(Int, 2) .* column
#
#     bd_pair = [orange_birth_r1, orange_death_r1]
#     repetitions = total_cols-(2*(column-1))
#     for m = 1:(repetitions)
#         plot!(x_vals, bd_pair; orange_kwargs...)
#     end
#
#     bd_pair = [blue_birth_r1, blue_death_r1]
#     repetitions2 = (2*(column))-total_cols
#     for m = 1:(repetitions2)
#         plot!(x_vals, bd_pair; blue_kwargs...)
#     end
#
#     for k=0:4
#         bd_pair = [orange_birth_r2, orange_death_r2] .+ k/10
#         repetitions3 = (column<=half_total_cols) ? (half_total_cols-k) : abs(column-total_cols)
#         for m = 1:(repetitions3)
#             plot!(x_vals, bd_pair; orange_kwargs...)
#         end
#
#         bd_pair = [blue_birth_r2, blue_death_r2] .+ k/10
#         repetitions4 = (column<=half_total_cols) ? column : (half_total_cols-k)
#         for m = 1:(repetitions4)
#             plot!(x_vals, bd_pair; blue_kwargs...)
#         end
#     end
#
# end
# display(barcodes_legend_canvas)
