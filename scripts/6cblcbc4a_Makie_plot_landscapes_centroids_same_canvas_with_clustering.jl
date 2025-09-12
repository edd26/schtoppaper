
#=
Creates a plot with coloured cluster, brain regions, vertical barcodes and landcsapes with their shuffling
using Makie library. Two possible backend could be used, by setting `cairo_backend` variable:
- if true, then CairoMakie us used
- else GLMakie is used
=#

import DrWatson: @quickactivate, srcdir, scriptsdir
@quickactivate "schtoppaper"
## ===-===-
"6cblcb_get_subjects_per_cluster_distributions.jl" |> scriptsdir |> include

"ClusterMergePlotting.jl" |> srcdir |> include
"ClusterMergePlottingMakie.jl" |> srcdir |> include
"DistributionsUtils.jl" |> srcdir |> include
"Centroids.jl" |> srcdir |> include

## ===-===-===-===-===-
using PersistenceLandscapes
#
# For dendrogram:
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
script_prefix = "6cblcbc4a"
@info "$(script_prefix): \tPlot shuffled landscapes's centroids in same canvas..."

if selected_dim == 0
    y_max = 0.1
    x_max = 0.2
elseif selected_dim == 1
    y_max = 0.06
    x_max = 0.3
elseif selected_dim == 2
    y_max = 0.05
    x_max = 1.0
else
    y_max = 0.025
    x_max = 1.0
end


## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
critical_values = get_BH_critical_values_form_landscapes(
    landscapes_distance,
    distances_distribution,
)
relevant_cluster_keys = critical_values[critical_values.AcceptedValue, "ClusterKeys"]
selected_clusters = parse.(Int, [split(k, "clust")[2] for k in relevant_cluster_keys])

all_keys = ["clust$(c|> get_cluster_split_name )" for (k, c) in hclust_m.clusters]
not_difference_cluster_keys = [k for k in all_keys if !in(k, relevant_cluster_keys)]
#
## ===-===-===-===-===-===-
# Colours setting
max_len = length(all_keys)
colours_dict = cgrad(:picasso, categorical=true, max_len, rev=false)
colours_dict_pt1 = cgrad(:cmyk, categorical=true, max_len, rev=false)
c1 = colours_dict_pt1[5:5:end]
colours_dict = cgrad([
        c1[1:3]...,
        :orange,
        c1[6:8]...,
        RGBf(([138, 186, 100] ./ 255)...),
        c1[10:end]...],
    categorical=true, max_len, rev=true)

groups_colouring = [k => PolyElement(color=color, strokecolor=:transparent)
                    for (k, color) in colors_subjects_with_names
] |> OrderedDict

legend_keys = ["hc" => "Healthy", "sch" => "Schizophrenia"] |> OrderedDict

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
total_histograms = distances_distribution |> length
total_cols = 1
total_rows = 1

points_per_cm = 28.3465 * 1.5
plt_width = 19.03 * points_per_cm
plt_height = 21.0 * points_per_cm 


final_plt_height = plt_height
final_plt_width = plt_width

# ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
f = Figure(size=(final_plt_width, final_plt_height));
fgl = GridLayout(f[1, 1]) # whole space

fgl_pland = GridLayout(fgl[1, 1]) # landscapes
fgl_clustering = GridLayout(fgl[2, 1]) # clustering

fga = GridLayout(fgl_pland[1, 1]) # plots

fgl_polar = CairoMakie.GridLayout(fgl_clustering[1, 1])
fgl_hmap = CairoMakie.GridLayout(fgl_clustering[1, 2])

# ===-===-
CairoMakie.rowsize!(fgl, 2, Relative(0.4))

# ===-===-===-===-===-===-===-===-===-===-
ax_hc = Axis(
    fga[1, 1],
    backgroundcolor=RGB(1.0, 0.8, 0.8),
    title="HC sublandscapes,\nsignificant",
    ylabel="Landscape\nheight [a. u.]",
    xtrimspine=(true, true),
)

ax_cent = Axis(
    fga[2, 1],
    backgroundcolor=RGB(1.0, 0.8, 0.8),
    title="Centroids in landscapes space\nsignificant",
    ylabel="Landscape\nheight [a. u.]",
    xtrimspine=(true, true),
)

ax_sch = Axis(
    fga[3, 1],
    backgroundcolor=RGB(1.0, 0.8, 0.8),
    title="SCH sublandscapes,\nsignificant",
    ylabel="Landscape\nheight [a. u.]",
    xlabel="Normalised filtration step",
    xtrimspine=(true, true),
)

# ===-
# not different
ax_hc_notsig = Axis(
    fga[1, 2],
    title="HC sublandscapes,\nno differences",
    xtrimspine=(true, true),
)
ax_cent_norsig = Axis(
    fga[2, 2],
    backgroundcolor=RGB(1.0, 1.0, 1.0),
    title="Centroids in landscapes space\nno difference",
    xlabel="Normalised filtration step",
    xtrimspine=(true, true),
)
ax_sch_notsig = Axis(
    fga[3, 2],
    title="SCH landcsapes,\n no differences",
    xlabel="Normalised filtration step",
    xtrimspine=(true, true),
)

hidespines!(ax_hc, :t, :r)
hidespines!(ax_sch, :t, :r)
hidespines!(ax_hc_notsig, :t, :r)
hidespines!(ax_sch_notsig, :t, :r)
hidespines!(ax_cent, :t, :r)
hidespines!(ax_cent_norsig, :t, :r)

ax_dict = Dict(
    "hc" => ax_hc,
    "sch" => ax_sch,
)

ax_dict_no_difference = Dict(
    "hc" => ax_hc_notsig,
    "sch" => ax_sch_notsig,
)
# ===-
# Landscvapes and histograms

ordered_permutation = #[mean(c.range) for (k, c) in hclust_m.clusters] |> sortperm
    [parse(Int, split(c, "clust")[2]) for c in relevant_cluster_keys] |> sortperm

x_hats = Dict([k => Float64[] for k in data_keys]...)
y_hats = Dict([k => Float64[] for k in data_keys]...)

x_hats2 = Dict([k => Float64[] for k in data_keys]...)
y_hats2 = Dict([k => Float64[] for k in data_keys]...)

################### >>>
cycle_keys_vectors = [relevant_cluster_keys, not_difference_cluster_keys,]

# Plot plandscapes
for (pland_axis_selection, cycles_selection) in zip([ax_dict, ax_dict_no_difference], cycle_keys_vectors)

    local_ordered_permutation = [parse(Int, split(c, "clust")[2]) for c in cycles_selection] |> sortperm

    for (cl_id, cluster_key) in cycles_selection[local_ordered_permutation] |> enumerate
        # ===-===-
        # set up row and col for plt
        # get landscape plots
        cluster_number = parse(Int, split(cluster_key, "clust")[2])
        @info "cluster_number $(cluster_number)"

        for (id, selected_key) in data_keys |> enumerate
            @info "\tselected_key $(selected_key)"


            pland = average_landscape_per_cluster[cluster_key][selected_key]
            x_vals, y_vals, x_hat, y_hat = get_landscape_centroid(pland)

            push!(x_hats[selected_key], x_hat)
            push!(y_hats[selected_key], y_hat)

            colours_group = colours_dict

            lines_kwargs = (color=colours_group[cluster_number], linewidth=3, alpha=0.7)
            CairoMakie.lines!(pland_axis_selection[selected_key], x_vals, y_vals; lines_kwargs...)
        end
    end
end

# Plot arrows
for (ax, clusters_selection) in zip([ax_cent, ax_cent_norsig], cycle_keys_vectors)
    @info ax

    local_ordered_permutation = 
        [parse(Int, split(c, "clust")[2]) for c in clusters_selection] |> sortperm
    for (cl_id, cluster_key) in clusters_selection[local_ordered_permutation] |> enumerate

        # ===-===-
        # set up row and col for plt
        # get landscape plots
        cluster_number = parse(Int, split(cluster_key, "clust")[2])
        @info "cluster_number $(cluster_number)"

        for (id, selected_key) in data_keys |> enumerate
            @info "\t selected_key $(selected_key)"

            pland = average_landscape_per_cluster[cluster_key][selected_key]
            x_vals, y_vals, x_hat, y_hat = get_landscape_centroid(pland)

            push!(x_hats[selected_key], x_hat)
            push!(y_hats[selected_key], y_hat)

            colours_group = colours_dict

            lines_kwargs = (color=colours_group[cluster_number], linewidth=3, alpha=0.6)

            # plot landsacpes
            CairoMakie.scatter!(ax, x_hat, y_hat, color=colors_subjects_with_names[selected_key])
        end
        arrow_head_x = [x_hats["sch"][end] - x_hats["hc"][end]]
        arrow_head_y = [y_hats["sch"][end] - y_hats["hc"][end]]
        arrow_tail_x = [x_hats["hc"][end]]
        arrow_tail_y = [y_hats["hc"][end]]


        selected_colour = nothing
        if clusters_selection == relevant_cluster_keys
            selected_colour = if cut_cluster[cl_id] == 1
                :green
            else
                :red
            end
        else
            selected_colour =
                :black
        end

        arrows!(
            ax,
            arrow_tail_x,
            arrow_tail_y,
            arrow_head_x,
            arrow_head_y,
            linewidth=5,
            lengthscale=0.8,
            color=(selected_colour, 0.8),
            alpha=0.8
        )
        f
    end
end

for ax in [ax_cent, ax_cent_norsig]
    low_x = -0.01
    high_x = 0.38
    x_ticks = 0:0.05:0.35
    makie_backend.xlims!(ax, low=low_x, high=high_x)
    makie_backend.ylims!(ax, low=0, high=0.065)
    ax.xticks = x_ticks
end

f
for ax in [ax_hc, ax_sch, ax_hc_notsig, ax_sch_notsig]
    low_x = -0.01
    high_x = 0.38
    x_ticks = 0:0.05:0.35

    makie_backend.xlims!(ax, low=low_x, high=high_x)
    makie_backend.ylims!(ax, low=0, high=0.065)

    ax.xticks = x_ticks
end



for ax in [ax_hc, ax_hc_notsig, ax_cent, ax_cent_norsig]
    hidexdecorations!(ax,
        ticks=false,
        grid=false,
        minorgrid=false,
        minorticks=false
    )
end
for ax in [ax_hc_notsig, ax_sch_notsig, ax_cent_norsig]
    hideydecorations!(ax,
        ticks=false,
        grid=false,
        minorgrid=false,
        minorticks=false
    )
end
f

angles_legend = [
    "Cluster 1",
    "Cluster 2",
]
fgl_legends2 = GridLayout(fga[:, end+1]) # legend

makie_backend.Colorbar(
    fgl_legends2[1, 1],
    limits=(1, max_len),
    size=10,
    colormap=colours_dict,
    label="Cycle index",
    tellheight=false,
    tellwidth=false,
    ticks=5:(ceil(Int, 0.1max_len)):max_len,
    flip_vertical_label=true
)

colsize!(fga, 3, Relative(0.01))
colgap!(fga, 2, 0)

f
########################3
# Second Figure
# ===-===-===-
ax_dendrogram = CairoMakie.Axis(fgl_hmap[1, 1],)
ax_ordered_hmap = CairoMakie.Axis(
    fgl_hmap[2, 1],
    aspect=AxisAspect(1),
    ylabel="Cycle index"
)
ax = PolarAxis(
    fgl_polar[1, 1],
    thetalimits=(
        -(4pi / 4),
        (4pi / 4),
    ),
    rtickangle=pi
)

# ===-===-===-
hm = CairoMakie.heatmap!(
    ax_ordered_hmap,
    vec_cosine_distance[ord, ord],
    colormap=:terrain,
    colorrange=(0, 2)
)
ax_ordered_hmap.yticks = ((1:total_landscapes), [split(l, "Cycle ")[2] for l in sorted_landsacpes_labels])

dendrogram_coloured!(ax_dendrogram, cos_distance_clustering;)
ax_dendrogram.ylabel = "Cluster\nheight"
CairoMakie.ylims!(ax_dendrogram, low=0)
ax_dendrogram.yticks = 0.5:0.5:1

# ===-===-
for (i, (r, theta)) in zip(all_r, all_theta) |> enumerate
    selected_colour = if cut_cluster[i] == 1
        :green
    else
        :red
    end
    total_points = length(all_r)
    lines!(
        ax,
        [theta for k in 1:total_points],
        [k for k in range(0, r, length=total_points)],
        linewidth=4,
        color=selected_colour,
        alpha=0.4,
    )
end

ax.rticklabelrotation = pi / 4

angles_legend = [
    "Cluster 1"
    "Cluster 2"
]
# ===-===-===-
fgl_legends1 = GridLayout(fgl_polar[end+1, :]) # legend
makie_backend.Legend(
    fgl_legends1[1, 1],
    [arrows_colouring],
    [angles_legend],
    ["Centroid vector angle"],
    nbanks=2,
    framevisible=false,
    tellheight=false,
    tellwidth=false,
    titlefont=:regular,
)


fgl_legends3 = GridLayout(fgl_hmap[2, 2]) # legend
makie_backend.Colorbar(
    fgl_legends3[1, 1],
    hm,
    label="Cosine distance",
    flip_vertical_label=true,
)

rowsize!(fgl_polar, 2, Relative(0.05))
colsize!(fgl_clustering, 1, Relative(0.53))

# join clistering and hmap
rowgap!(fgl_hmap, 1, 2)
rowsize!(fgl_hmap, 1, Relative(0.2))

linkxaxes!(ax_dendrogram, ax_ordered_hmap)
hidespines!(ax_dendrogram, :t, :r)
hidexdecorations!(ax_dendrogram)
hidexdecorations!(ax_ordered_hmap)


for (i, fgl_local) in enumerate([fgl_pland, fgl_polar, fgl_hmap])

    selected_padding =
        if i == 1
            (0, 70, 10, 0)
        else
            (0, 70, -15, 0)
        end

    label = "$(('a':'z')[i]))"
    Label(fgl_local[1, 1, TopLeft()], label,
        fontsize=18,
        padding=selected_padding,
        halign=:right,
        tellwidth=false
    )
end


f

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
@info "$(script_prefix): Saving centroids direction for all keys..."

## ===-===-===-
ENV["GKSwstype"] = "100"

## ===-===-===-
plots_6cblcbc3_dir(args...) = plotsdir("section6", script_prefix * "_landscapes_centroid_single_canvas", args...)
landscapes_prefix = ""

## ===-===-===-
dim = selected_dim
p = p_value
pathargs = (1,)
if manual_clustering > 0
    splits = manual_clustering
    names_to_save = @dict savename joined_keys total_clusters linkage dim p total_shuffles splits minimal_height normalisation_type
    savename_val = savename(names_to_save, allowedtypes=(DrWatson.default_allowed(names_to_save)..., NormalisationType))
    pathargs = (
        "n-times_split",
        "joined_keys=$(joined_keys)",
        "dim$(selected_dim)",
        "linkage=$(linkage)",
        "total_clust=$(total_clusters)_total_shuffles=$(total_shuffles)",
        "normalisation_type=$(normalisation_type)",
    )
else
    names_to_save = @dict joined_keys cluster_height min_clusters total_clusters linkage dim p total_shuffles normalisation_type
    savename_val = savename(names_to_save, allowedtypes=(DrWatson.default_allowed(names_to_save)..., NormalisationType))
    pathargs = (
        "joined_keys=$(joined_keys)",
        "dim$(selected_dim)",
        "linkage=$(linkage)",
        "total_clust=$(total_clusters)",
        "total_shuffles=$(total_shuffles)",
        "normalisation_type=$(normalisation_type)",
    )
end


if (("HCP_test_orig" in data_keys) || ("dHCP_test_orig" in data_keys)) && do_HCP_trim
    savename_val *= "_HCP_trim=$(do_HCP_trim)_ord_limit=$(ordering_limit)"
end

if do_vector_extension
    savename_val *= "_vector_extension=$(do_vector_extension)"
end

if limit_popularity != 0
    pathargs = (pathargs..., "popularity_thr=$(limit_popularity)")
    savename_val *= "_popularity_thr=$(limit_popularity)"
end

if repeat_cycles
    pathargs = (pathargs..., "repeat_cycles")
    savename_val *= "_repeat_cycles=$(repeat_cycles)"
end

## ===-===-===-
landscape_plot_name = plots_6cblcbc3_dir(pathargs..., script_prefix * "$(landscapes_prefix)_$(savename_val).png")
safesave(landscape_plot_name, f)

landscape_plot_name_pdf = plots_6cblcbc3_dir(pathargs..., "pdf", script_prefix * "$(landscapes_prefix)_$(savename_val).pdf")
safesave(landscape_plot_name_pdf, f, pdf_version="1.5")
save(landscape_plot_name_pdf, f, pdf_version="1.5")

@info "$(script_prefix): Saved as" landscape_plot_name

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
do_nothing = "ok"
