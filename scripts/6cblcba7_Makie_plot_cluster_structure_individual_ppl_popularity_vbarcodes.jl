
#=
Creates a plot with coloured cluster, brain regions, vertical barcodes and landcsapes with their shuffling
using Makie library. Two possible backend could be used, by setting `cairo_backend` variable:
- if true, then CairoMakie us used
- else GLMakie is used
=#

import DrWatson: @quickactivate, srcdir, scriptsdir
@quickactivate "schtoppaper"
## ===-===-
"6cblcb_2_get_subjects_per_cluster_distributions_simplified.jl" |> scriptsdir |> include

"ClusterMergePlotting.jl" |> srcdir |> include
"ClusterMergePlottingMakie.jl" |> srcdir |> include
"DistributionsUtils.jl" |> srcdir |> include

## ===-===-===-===-===-
using PersistenceLandscapes
using StatsBase

## ===-
cairo_backend = true
DO_TWO_TAILED = false
NTH_PERCENTILE = 0.95
TOTAL_BINS = 100

if manual_clustering > 0
    total_clusters = 2^manual_clustering
else
    total_clusters = clusters_keys |> length
end

if cairo_backend
    using CairoMakie
    makie_backend = CairoMakie
else
    using GLMakie
    makie_backend = GLMakie
end
## ===-===-===-===-===-
script_prefix = "6cblcba7"
@info "$(script_prefix): \tPlot shuffled landscapes's distributions ..."

if selected_dim == 0
    y_max = 0.1
    x_max = 0.2
elseif selected_dim == 1
    y_max = 0.08
    x_max = 1.0
elseif selected_dim == 2
    y_max = 0.05
    x_max = 1.0
else
    y_max = 0.025
    x_max = 1.0
end

## ===-===-===-===-===-===-
# Colours setting
rng = MersenneTwister(1234);
total_cluster_name = get_max_cluster_index(hclust_m)
col_palette = palette(:grays, total_cluster_name; rev=true);
shuffled_colours = col_palette[shuffle(rng, 1:total_cluster_name)];

colors_brain_regions = [RGBf(k) for k in colours_palete];
colors_brain_regions[4] =
    RGB(
        colors_brain_regions[4].r * 1.1,
        colors_brain_regions[4].g * 0.9,
        colors_brain_regions[4].b * 0.1,
    )

if "HCP_1" in data_keys && "hc" in data_keys && "sch" in data_keys
    push!(colors_subjects, data_colours_bank["COBRE"])
end


brain_colouring = [
    PolyElement(color=color, strokecolor=:transparent)
    for color in colors_brain_regions
];


## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
hc_popularity = sum(subject_presence_in_cycle[1:44, :] .|> !ismissing, dims=1)
sch_popularity = (-1) .* sum(subject_presence_in_cycle[45:end, :] .|> !ismissing, dims=1)
total_cycles = size(subject_presence_in_cycle, 2)
tick_step = 0
if total_cycles > 5000
    tick_step = 500
elseif total_cycles > 2500
    tick_step = 200
else
    tick_step = 50
end


raster_plt_ylabel = ""
popularity_plt_ylabel = ""
popularity2_plt_ylablel = ""
raster_plt_colopmap = []

pop_base_y_label = "Popularity"
if "hc" in data_keys && "sch" in data_keys
    raster_plt_ylabel = "AAL2 atlas index"
    popularity_plt_ylabel = "Participants'\nindex"
    if "HCP_1" in data_keys
        popularity2_plt_ylablel = "$(pop_base_y_label)\n   HCP ->\n<- COBRE"
    else
        popularity2_plt_ylablel = "$(pop_base_y_label)\n   SCH ->\n<- HC"
    end

    raster_plt_colopmap = colors_brain_regions
else
    raster_plt_ylabel = "Model nodes"
    popularity_plt_ylabel = "Sample index"
    popularity2_plt_ylablel = "$(pop_base_y_label)\n GROUP 1 ->\n<- GROUP 2"
    raster_plt_colopmap = cgrad(:greys, 8, categorical=true, rev=true)[2:end]
end

# Source: 
# https://docs.makie.org/stable/explanations/figure#Figure-size-and-units
# 23.12cm
# 15.66cm
# 1cm = 28.3465
# 28.3465 * 23.12
# 28.3465 * 15.66

all_size_set = ["no-legend", "small"]
size_set = ["no-legend", "small"][1]
for size_set in all_size_set
    points_per_cm = 0
    points_per_cm = 28.3465 * 2
    plt_width = floor(Int, 18 * points_per_cm)
    plt_height = floor(Int, 13 * points_per_cm)

    if size_set == "no-legend"
        plt_height -= 30
    end

    f = Figure(size=(plt_width, plt_height,), pt_per_unit=1)

    # Setting up GridLayouts
    fgl = f[1, 1] = GridLayout() # whole space

    ga = fgl[1, 1] = GridLayout() # main plots
    ga1 = ga[1, 1] = GridLayout() # dendrogram
    ga2 = ga[2, 1] = GridLayout() # brain regions
    ga3 = ga[3, 1] = GridLayout() # popularity
    ga4 = ga[4, 1] = GridLayout() # popularity
    ga5 = ga[5, 1] = GridLayout() # vertical barcodes

    gb = Nothing
    gb2 = Nothing
    gb3 = Nothing
    if size_set == "no-legend"
        print("No-legend option")
    else
        gb = fgl[2, 1] = GridLayout() # legends
        gb2 = gb[1, 1] = GridLayout() # brain regions legend
        gb3 = gb[1, 2] = GridLayout() # subjectrs legend
    end

    gc = fgl[1, 0] = GridLayout() # landscapes and histograms
    gc_grids = [gc[k, 1] |> GridLayout for k in 1:5]
    subsubscripts = ["i", "ii", "iii", "iv", "v"]

    for (ref, grid,) in zip(subsubscripts, gc_grids)
        label = "$(ref))"
        CairoMakie.Label(
            grid[1, 1],
            label,
            tellheight=false,
            valign=:top,
            fontsize=16
        )
    end

    rowgap!(fgl, 0)

    ## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
    # Panel A - cluster and raster
    # ===-
    # Plot dendrogram
    ax_dendrogram = Axis(ga1[1, 1])

    dendrogram_coloured!(ax_dendrogram, hclust_m; col_palette=shuffled_colours)

    ax_dendrogram.xticklabelrotation = pi / 4

    max_height = max(selected_cluster.height...)
    makie_backend.ylims!(ax_dendrogram, 0, 1.1max_height)
    makie_backend.xlims!(ax_dendrogram, 0, length(merged_cluster))


    ax_dendrogram.yticks = 0:100:1.1max_height
    ax_dendrogram.ylabel = "Cluster height\n[a.u.]"

    hidexdecorations!(ax_dendrogram,)
    hideydecorations!(ax_dendrogram, label=false, grid=false)
    f

    ## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
    # Preparing brain regions plots

    cycles_presence = large_scale_network_matrix_preprocessing(regions_in_cycles |> replace_zeros_with_missing,
        yeo7_networks, brain_regions_info
    )
    regions_presence = cycles_presence[new_regions_order, selected_cluster.order]' |> Matrix

    values_x = Float64[]
    labels_x = String[]

    for (n, c) in hclust_m.clusters
        clust_name = c |> get_cluster_split_name

        left = c.range[1]
        right = c.range[end]
        push!(values_x, left + 0.5)
        push!(labels_x, "$(clust_name)")
    end

    ## ===-===-===-===-===-===-===-===-
    # plotting brian regions
    ax_brain_regions = Axis(ga2[1, 1])

    # if not hc or sch, elemnts are not regions
    if "hc" in data_keys && "sch" in data_keys
        nothing
    else
        not_missing = findall(x -> !ismissing(x), regions_presence)
        missing_elements = findall(x -> ismissing(x), regions_presence)
        regions_presence[not_missing] .= 1
    end
    regions_presence = regions_presence |> Matrix{Union{Float32,Missing}}

    makie_backend.heatmap!(
        ax_brain_regions,
        regions_presence,
        colorrange=(1, 8),
        alpha=0.9,
        colormap=raster_plt_colopmap,
    )

    # ===-===-===-
    # Add Axis 
    ax_brain_regions.ylabel = raster_plt_ylabel
    hideydecorations!(ax_brain_regions, grid=true, label=false)

    ax_brain_regions.xticks = 0:tick_step:total_cycles
    hidexdecorations!(ax_brain_regions)#, grid=false)

    ## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
    # Popularity plot
    ax_popularity = Axis(ga3[1, 1])

    makie_backend.heatmap!(
        ax_popularity,
        subject_presence_in_cycle[:, selected_cluster.order]',
        colormap=colors_subjects[1:length(data_keys)]
    )

    ppl_step = 5
    ax_popularity.ylabel = "Subject index"
    if length(data_keys) == 1
        vals_range = 1:ppl_step:total_matrices
        values_py = [v for v in vals_range]
        labels_py = [
            ["$(data_keys[1])_$(k)" for k in 1:total_matrices]...,
        ][vals_range]
    else
        vals_range = 1:ppl_step:2total_matrices
        values_py = [v for v in vals_range]
        labels_py = [
            ["$(data_keys[1])_$(k)" for k in 1:total_matrices]...,
            ["$(data_keys[2])_$(k)" for k in 1:total_matrices]...,
        ][vals_range]
    end

    ax_popularity.yticks = (values_py, labels_py)
    ax_popularity.ylabel = popularity_plt_ylabel
    ax_popularity.xticks = 0:tick_step:total_cycles
    hidexdecorations!(ax_popularity)
    hideydecorations!(ax_popularity, label=false)

    ## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
    # Popularity plot
    ax_popularity2 = Axis(ga4[1, 1])


    if length(data_keys) > 2 && "hc" in data_keys && "sch" in data_keys
        sch_popularity = sum(subject_presence_in_cycle[45:88, :] .|> !ismissing, dims=1)
        hc_popularity = sum(subject_presence_in_cycle[1:44, :] .|> !ismissing, dims=1)

        hcp_popularity = sum(subject_presence_in_cycle[89:end, :] .|> !ismissing, dims=1)

        plot_data1 = hcp_popularity
        plot_data2 = (-1) .* (sch_popularity + hc_popularity)
    else
        sch_popularity = sum(subject_presence_in_cycle[45:end, :] .|> !ismissing, dims=1)
        hc_popularity = (-1) .* sum(subject_presence_in_cycle[1:44, :] .|> !ismissing, dims=1)

        plot_data1 = sch_popularity
        plot_data2 = hc_popularity
    end

    if "HCP_1" in data_keys && "hc" in data_keys && "sch" in data_keys
        bar_top_colour = data_colours_bank["HCP_1"]
        bar_bot_colour = data_colours_bank["COBRE"]
    else
        bar_top_colour = colors_subjects[end]
        bar_bot_colour = colors_subjects[1]
    end
    barplot!(
        ax_popularity2,
        1:total_cycles,
        plot_data1[selected_cluster.order],
        color=bar_top_colour
    )
    barplot!(
        ax_popularity2,
        1:total_cycles,
        plot_data2[selected_cluster.order],
        color=bar_bot_colour
    )

    begin
        "Preparation for popularity analysis"
        reordered_sch_popularity = sch_popularity[selected_cluster.order]
        reordered_hc_popularity = hc_popularity[selected_cluster.order] .* (-1)

        total_ppl_per_cycle = reordered_sch_popularity .+ reordered_hc_popularity

        differences_in_popularity = abs.(reordered_sch_popularity .- reordered_hc_popularity)

        found_in_at_least_2 = vcat(findall(x -> x > 1, reordered_hc_popularity)...,
                                  findall(x -> x > 1, reordered_sch_popularity)...) |> unique |> sort

        differences_in_popularity = abs.(reordered_sch_popularity[found_in_at_least_2] .- reordered_hc_popularity[found_in_at_least_2])
    end


    # ===-
    # Axis management
    if all([d in data_keys for d in ["HCP_1", "hc", "sch"]])
        y_pop_max = 88
    else
        y_pop_max = 44
    end

    makie_backend.xlims!(
        ax_popularity2,
        low=0.0,
        high=size(ordered_bd_type_split_updated[data_keys[1]], 1)
    )
    makie_backend.ylims!(
        ax_popularity2,
        low=-y_pop_max,
        high=y_pop_max
    )

    ax_popularity2.ylabel = popularity2_plt_ylablel
    ytick_step = y_pop_max ÷ 2
    ax_popularity2.yticks = (-y_pop_max:ytick_step:y_pop_max, ["$(abs(k))" for k in (-y_pop_max:ytick_step:y_pop_max)])

    ax_popularity2.xticks = 0:tick_step:total_cycles
    ax_popularity2.xticklabelrotation = pi / 4
    hidexdecorations!(ax_popularity2, grid=true, label=true, ticks=true, ticklabels=true)
    ## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
    # Plot vertical barcodes
    ax_vert_barcodes = Axis(ga5[1, 1])
    mplot_vertical_barcodes!(
        ax_vert_barcodes,
        ordered_bd_type_split_updated,
        selected_cluster.order;
        alpha_lvl=0.45,
        color_bank=colors_subjects
    )

    # ===-
    # Axis management
    makie_backend.xlims!(
        ax_vert_barcodes,
        low=0.5,
        high=size(ordered_bd_type_split_updated[data_keys[1]], 1) + 0.5
    )

    if selected_dim == 0
        y_max = 0.05
    elseif selected_dim == 1
        y_max = 0.19
    elseif selected_dim == 1
        y_max = 0.5
    else
        y_max = 1.0
    end


    makie_backend.ylims!(
        ax_vert_barcodes,
        low=0,
        high=y_max
    )

    ax_vert_barcodes.ylabel = "Normalized\nlifetime"
    ax_vert_barcodes.xlabel = "Cycle index"
    ax_vert_barcodes.xticklabelrotation = pi / 4
    ax_vert_barcodes.xticks = 0:tick_step:total_cycles

    hidexdecorations!(ax_vert_barcodes, label=false, ticklabels=false)
    
    ## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
    # Add legends and colorbars

    yeo_full_names = [
        "Default" => "Default",
        "Limbic" => "Limbic",
        "DorsAttn" => "Dorsal Attention",
        "Vis" => "Visual",
        "Cont" => "Frontoparietal",
        "SalVentAttn" => "Ventral Attention",
        "Subcortical" => "Subcortical",
        "SomMot" => "Somatomotor",] |> Dict

    subject_keys = Dict(
        "hc" => "Healthy",
        "sch" => "Schizophrenia",
        "HCP_test_orig" => "Non inversed HCP",
        "HCP_test_" => "Inversed HCP",
        "HCP_1" => "HCP",
    )

    if size_set != "no-legend"
        if "hc" in data_keys && "sch" in data_keys
            makie_backend.Legend(
                gb2[1, 1],
                [brain_colouring],
                [[yeo_full_names[k] for (k, v) in ls_region_to_color]],
                ["Brain regions",],
                nbanks=2,
                orientation=:horizontal,
                framevisible=false,
            )
        end

        legend_keys = []
        for (i, k) in enumerate(data_keys)
            if haskey(subject_keys, k)
                push!(legend_keys, subject_keys[k])
            else
                push!(legend_keys, "Group $(i)")
            end
        end
        if "HCP_1" in data_keys && "hc" in data_keys && "sch" in data_keys
            push!(legend_keys, "COBRE")
        end

        groups_colouring = [PolyElement(color=color, strokecolor=:transparent)
                            for color in colors_subjects
        ]
        makie_backend.Legend(
            gb3[1, 1],
            [groups_colouring],
            [legend_keys],
            ["Groups",],
            orientation=:horizontal,
            framevisible=false,
            nbanks=2
        )
    end # if size_set 

    # ===-===-===-===-===-===-===-
    # set cluster and brain regions spacing

    size_ga1 = 0.15
    size_ga2 = 0.35
    size_ga3 = 0.20
    size_ga4 = 0.10
    size_ga5 = 0.20

    makie_backend.rowgap!(ga, 10)
    makie_backend.colgap!(ga1, 10)
    makie_backend.rowgap!(gc, 10)

    for gx in [ga, gc]
        makie_backend.rowsize!(gx, 1, Relative(size_ga1))
        makie_backend.rowsize!(gx, 2, Relative(size_ga2))
        makie_backend.rowsize!(gx, 3, Relative(size_ga3))
        makie_backend.rowsize!(gx, 4, Relative(size_ga4))
        makie_backend.rowsize!(gx, 5, Relative(size_ga5))
    end
    makie_backend.colgap!(fgl, 5)

    if size_set !== "no-legend"
        makie_backend.colsize!(gb, 1, Relative(0.6))
        makie_backend.rowsize!(fgl, 1, Relative(0.9))
    end
    f

    ## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
    @info "$(script_prefix): Saving persistence landscapes for each key..."

    ## ===-===-===-
    ENV["GKSwstype"] = "100"

    ## ===-===-===-
    script_subname = "dendrogram_structure_popularity"
    plots_6cblcbaa4_dir(args...) = plotsdir("section6", script_prefix * "_$(script_subname)", args...)
    landscapes_prefix = ""

    ## ===-===-===-
    dim = selected_dim
    p = p_value
    savename_val = @savename joined_keys cluster_height min_clusters linkage dim normalisation_type
    pathargs = (
        "joined_keys=$(joined_keys)",
        "dim$(selected_dim)",
        "linkage=$(linkage)",
        "$(size_set)",
        "normalisation_type=$(normalisation_type)"
    )

    if limit_popularity != 0
        pathargs = (pathargs..., "popularity_thr=$(limit_popularity)")
        savename_val *= "_popularity_thr=$(limit_popularity)"
    end

    if repeat_cycles
        pathargs = (pathargs..., "repeat_cycles")
        savename_val *= "_repeat_cycles=$(repeat_cycles)"
    end

    ## ===-===-===-
    landscape_plot_name = plots_6cblcbaa4_dir(pathargs..., script_prefix * "$(landscapes_prefix)_$(savename_val).png")
    safesave(landscape_plot_name, f)

    landscape_plot_name_pdf = plots_6cblcbaa4_dir(pathargs..., "pdf", script_prefix * "$(landscapes_prefix)_$(savename_val).pdf")
    safesave(landscape_plot_name_pdf, f)

    @info "$(script_prefix): Saved as" landscape_plot_name
end

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
do_nothing = "ok"
