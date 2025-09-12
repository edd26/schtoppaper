using DrWatson
@quickactivate "schtoppaper"

"1g_landcscapes_samples.jl" |> scriptsdir |> include

using CairoMakie
# ===-===-
x1 = 42
x2 = 43
x3 = 44
x4 = 45

mat1 = [
    0 1 10 16 1 5 7
    1 0 11 1 1 6 17
    10 11 0 12 13 19 18
    16 1 12 0 1 8 1
    1 1 13 1 0 9 15
    5 6 19 8 9 0 21
    7 17 18 1 15 21 0
]

mat1 = [
    0 1 x1 4 5 9 13 17
    1 0 2 x2 6 10 14 18
    x1 2 0 3 7 11 15 19
    4 x2 3 0 8 12 16 20
    5 6 7 8 0 21 x3 24
    9 10 11 12 21 0 22 x4
    13 14 15 16 x3 22 0 23
    17 18 19 20 24 x4 23 0
]
mat_1 = [
    0 1 42 4 9 13 21 28 31 34
    1 0 2 43 10 14 26 29 32 35
    42 2 0 3 11 15 27 30 33 36
    4 43 3 0 12 16 37 38 39 40
    9 10 11 12 0 41 22 23 24 25
    13 14 15 16 41 0 17 18 19 20
    21 26 27 37 22 17 0 5 44 8
    28 29 30 38 23 18 5 0 6 45
    31 32 33 39 24 19 44 6 0 7
    34 35 36 40 25 20 8 45 7 0
]
mat_2 = [
    0 1 42 4 9 13 25 28 31 34
    1 0 2 43 10 14 26 29 32 35
    42 2 0 3 11 15 27 30 33 36
    4 43 3 0 12 16 37 38 39 40
    9 10 11 12 0 41 21 22 23 24
    13 14 15 16 41 0 17 18 19 20
    25 26 27 37 21 17 0 5 44 8
    28 29 30 38 22 18 5 0 6 45
    31 32 33 39 23 19 44 6 0 7
    34 35 36 40 24 20 8 45 7 0
]


selected_mat = copy(mat_2)
total_rows, total_cols = size(selected_mat)

max_val = max(selected_mat...)

do_reverse = false
if do_reverse
    mat = selected_mat

    selected_mat = .-mat .+ max_val

    total_rows = size(mat, 1)
    diagonal_indices = [CartesianIndex(k, k) for k in 1:total_rows]
    selected_mat[diagonal_indices] .= 0

    last_step = max_val
    x_max_limit = max_val
else
    last_step = 38
    x_max_limit = 46 # last_step
end


# Run eirene computations
min_dim = 0
max_dim = 4
C = eirene(selected_mat, maxdim=max_dim, model="vr",)

data_mat = miss_elements_on_diagonal(selected_mat)


barcodes = get_barcodes(C, max_dim; min_dim=min_dim)
bettis = get_vectorized_bettis(C, max_dim; min_dim=min_dim)

topo_features = Dict()
topo_features["norm_barcodes"] = TopologyPreprocessing.get_normalised_barcodes(barcodes, bettis)
topo_features["barcodes"] = get_barcodes(C, max_dim; min_dim=min_dim)
topo_features["bettis"] = get_vectorized_bettis(C, max_dim; min_dim=min_dim)


function plt_barcodes!(ax, barcodes, dim_range=1:size(barcodes, 1); sort_by_birth=true, colours_palette=get_bettis_color_palete(min_dim=min_dim), last_step::Union{Int,Nothing}=nothing, max_x::Union{Float64,Nothing}=nothing)
    max_visible_step = max(max([[m for m in k if !isinf(m)] for k in barcodes]...)...)
    max_step =
        if isnothing(last_step)
            max_visible_step
        else
            last_step
        end

    last_y = 0
    for (b, bd_data) in barcodes |> enumerate
        # check if barcodes have non-zero length
        all_non_zero_barcodes = findall(x -> x != 0, bd_data[:, 2] - bd_data[:, 1])
        if all_non_zero_barcodes |> isempty
            continue
        else
            bd_data = bd_data[all_non_zero_barcodes, :]
        end

        current_dim = dim_range[b]
        all_infs = findall(x -> isinf(x), bd_data)
        if !isempty(all_infs)
            @warn "Infinity interval is replaced with a max step set to $(max_step)"
            bd_data[all_infs] .= max_step
        end

        total_lines = length(bd_data[:, 1])
        y_vals = collect(1:total_lines) .+ last_y
        last_y = max(y_vals...)

        if sort_by_birth
            if current_dim == 0
                barcodes_sorting = sortperm(bd_data[:, 2])
            else
                barcodes_sorting = sortperm(bd_data[:, 1])
            end
        else
            barcodes_sorting = 1:length(bd_data[:, 2])
        end

        lows = bd_data[barcodes_sorting, 1]
        highs = bd_data[barcodes_sorting, 2]
        rangebars!(ax, y_vals, lows, highs, color=colours_palette[b],
            whiskerwidth=10, direction=:x)
        # @info "Should print infinity at ($(max_visible_step), $(last_y)), current_dim=$(current_dim) is inf $(!isempty(all_infs))"
        if current_dim == 0 && !isempty(all_infs)
            # text!(ax, max_visible_step, last_y-1, text="Inf")#,align=(:right, :bottom))
            # text!(ax, 2, last_y - 1, text="Inf")#,align=(:right, :bottom))
            # @info "Printing infinity at ($(max_visible_step), $(last_y))"
            inf_marker_position =
                if isnothing(max_x)
                    max_visible_step
                else
                    max_x
                end
            # text!(ax, inf_marker_position, last_y, text=L"\infty", align=(:left, :bottom))
            f
        end
    end

end

## ===-===-
total_items = ceil(Int, (total_rows - 1) * total_rows / 2)
# colours_palette = cgrad(:Paired_12, total_items, categorical=true, rev=true);
# colours_palette = cgrad(:plasma, total_items, categorical=true, rev=false);
colours_palette = cgrad(:roma, total_items, categorical=true, rev=true);

angles = range(0, stop=2π, length=total_rows + 1)
x_coords = sin.(angles)
y_coords = cos.(angles)

# special shape
#           A .    B .    C .   D .     E .    F .     G .     H .     I .     J
x_coords = [947.7, 567.8, 93.0, 464.6, 1230.6, 1230.6, 1516.8, 1888.4, 2371.5, 1991.7]
y_coords = [749.9, 858.3, 803.5, 709.8, 1436.1, 114.4, 763.3, 669.5, 709.6, 818]
# x_coords = [947.7, 569.1, 93.0, 465.8, 1230.6, 1230.6, 1518.8, 1888.4, 2371.5, 1991.7]
# y_coords = [749.9, 967.3, 891.5, 818.6, 1436.1, 114.4, 914.7, 669.5, 709.6, 818]
#           A .  B .  C .  D .  E .  F .  G .  H .  I    J
# x_coords = [420, 200, 35, 190, 500, 500, 810, 580, 800, 965]
# y_coords = [360, 620, 320, 140, 750, 50, 140, 360, 620, 320]

alphabet = collect('A':'Z')[1:total_rows]
xticks = yticks = string.(collect(alphabet))

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
# 1cm <=> 28.3465 pt
# 6.3 * 28.3465
points_per_cm = 28.3465

plt_width = 17.7 * points_per_cm * 1.5
plt_height = 20.00 * points_per_cm * 1.5
f = CairoMakie.Figure(size=(plt_width, plt_height,), pt_per_unit=1);
f_main = GridLayout(f[1, 1])
fgl = GridLayout(f_main[1, 1])

# ===-
# Plot landscapes
gb = GridLayout(fgl[1, 2])
gc = GridLayout(fgl[2, 2])
ga = GridLayout(fgl[1:2, 1])
gd = GridLayout(fgl[3, :])

ax_land = CairoMakie.Axis(gb[1, 1], xminorticks=IntervalsBetween(2),)
ax_barcodes = CairoMakie.Axis(gc[1, 1], xminorticks=IntervalsBetween(2))
ax_hmap = CairoMakie.Axis(
    ga[1, 1],
    xticks=(1:total_rows, xticks),
    yticks=(1:total_rows, yticks),
    aspect=AxisAspect(1),
    title = "Connectivity matrix",
)

# ===-===-
# Create persistence landscapes

dim_x_landscape = PersistenceLandscape[]
for selected_dim = min_dim:max_dim
    sdim_index = (min_dim == 0) ? (selected_dim + 1) : (selected_dim)

    land_data = topo_features["barcodes"][sdim_index]
    total_bars = size(land_data, 1)

    push!(dim_x_landscape,
        [MyPair(land_data[k, :]...) for k in 1:total_bars] |> PersistenceBarcodes |> PersistenceLandscape
    )
end


colors = get_bettis_color_palete(min_dim=min_dim);
for p = 1:(max_dim+1)
    pl1 = dim_x_landscape[p]
    max_layers = size(pl1.land, 1)

    for k = 1:max_layers
        peaks_position, peaks = PersistenceLandscapes.get_peaks_and_positions(pl1.land[k])

        CairoMakie.lines!(
            ax_land,
            peaks_position,
            peaks;
            color=(colors[p], 1 / k),
            linewidth=6
        )
    end
end

# ===-
# Plot barcodes
plt_barcodes!(ax_barcodes, topo_features["barcodes"], min_dim:max_dim;
    last_step=45,
    max_x=43.5)

# ===-
# Plot Hmap

hmap = CairoMakie.heatmap!(ax_hmap, data_mat, colormap=colours_palette,)
for i in 1:total_rows, j in i:total_rows
    if i == j
        continue
    end
    txtcolor = (9 > data_mat[i, j]) || (data_mat[i, j] >= 37) ? :white : :black

    text!(ax_hmap, "$(data_mat[i,j])", position=(i, j), color=txtcolor, align=(:center, :center))
    if i != j
        text!(ax_hmap, "$(data_mat[j,i])", position=(j, i), color=txtcolor, align=(:center, :center))
    end
end
ax_hmap.xticklabelalign = (:right, :center)

# ===-===-===-===-===-===-===-===-===-===-===-===-
# Scaffolds
last_row = 0
col_index = 0
images_per_row = 8

function add_annotations!(ax, x_coords, y_coords)
    coords_range = [3, 4, 7, 8, 9]
    annotations!(
        ax,
        xticks[coords_range],
        [Point(x_coords[n], y_coords[n]) for n in coords_range],
        align=(:right, :top),
        fontsize=11,
        strokecolor=:white,
    )
    coords_range2 = [1, 2, 5, 6, 10]
    annotations!(
        ax,
        xticks[coords_range2],
        [Point(x_coords[n], y_coords[n]) for n in coords_range2],
        align=(:left, :bottom),
        fontsize=11,
        strokecolor=:white,
    )

    CairoMakie.xlims!(ax, low=-400, high=1.2max(x_coords...))
    CairoMakie.ylims!(ax, low=-10, high=1.2max(y_coords...))
end
function add_dots!(gd, last_row, col_index; images_per_row=8)
    if col_index == images_per_row # && col_index!=1
        last_row += 1
        col_index = 1
    else
        col_index += 1
    end
    CairoMakie.Label(gd[last_row, col_index], "...", tellwidth=false, tellheight=false)

    return last_row, col_index
end

function add_structure_at_position(gd, last_row, col_index, k, selected_mat; images_per_row=8)
    if col_index == images_per_row
        last_row += 1
        col_index = 1
    else
        col_index += 1
    end

    scaffold_matrix = copy(selected_mat)
    scaffold_matrix[findall(x -> x > k, selected_mat)] .= 0

    ax = CairoMakie.Axis(gd[last_row, col_index], title="$(k)")
    CairoMakie.scatter!(ax, x_coords, y_coords)
    append_scaffold(
        ax,
        scaffold_matrix,
        x_coords,
        y_coords,
        colours_palette;
        data_matrix=replace(data_mat, missing => 0,),
        scaffold_kwargs...
    )
    add_annotations!(ax, x_coords, y_coords)
    hidedecorations!(ax,)
    hidespines!(ax)
    return last_row, col_index
end

scaffold_kwargs = (lw=3, alpha=0.6)

for k in 1:23 # to accomodate dimension 4 structure
    if k % images_per_row == 1# && col_index!=1
        global last_row += 1
        global col_index = 1
    else
        global col_index += 1
    end

    scaffold_matrix = copy(selected_mat)
    scaffold_matrix[findall(x -> x > k, selected_mat)] .= 0

    ax = CairoMakie.Axis(gd[last_row, col_index], title="$(k)")
    CairoMakie.scatter!(ax, x_coords, y_coords)


    append_scaffold(
        ax,
        scaffold_matrix,
        x_coords,
        y_coords,
        colours_palette;
        data_matrix=replace(data_mat, missing => 0,),
        scaffold_kwargs...
    )

    add_annotations!(ax, x_coords, y_coords)

    hidedecorations!(ax,)
    hidespines!(ax)
end

k = 24
last_row, col_index = add_structure_at_position(gd, last_row, col_index, k, selected_mat)
k = 25
last_row, col_index = add_structure_at_position(gd, last_row, col_index, k, selected_mat)


last_row, col_index = add_dots!(gd, last_row, col_index; images_per_row=images_per_row)
k = 28
last_row, col_index = add_structure_at_position(gd, last_row, col_index, k, selected_mat)

last_row, col_index = add_dots!(gd, last_row, col_index; images_per_row=images_per_row)
k = 40
last_row, col_index = add_structure_at_position(gd, last_row, col_index, k, selected_mat)
k = 41
last_row, col_index = add_structure_at_position(gd, last_row, col_index, k, selected_mat)

last_row, col_index = add_dots!(gd, last_row, col_index; images_per_row=images_per_row)
k = 45
last_row, col_index = add_structure_at_position(gd, last_row, col_index, k, selected_mat)

f

# ===-===-===-===-
# Add plot annotations

for (label, layout) in zip(["a)", "b)", "c)", "d)"], [ga, gb, gc, gd])
    Label(layout[1, 1, TopLeft()], label,
        fontsize=18,
        padding=(0, 2, 2, 0),
        halign=:right)
end

## ===============
# General cleanup

ax_hmap.yreversed = true

hidexdecorations!(ax_hmap, ticks=false, ticklabels=false)
hideydecorations!(ax_hmap, ticks=false, ticklabels=false)
hidexdecorations!(ax_land, grid=false)
hideydecorations!(ax_land, grid=true, label=false)
hideydecorations!(ax_barcodes, grid=true, label=false)
hidespines!(ax_land, :t, :r, :l)
hidespines!(ax_barcodes, :t, :r, :l)
rowgap!(fgl, 10)

# set sections spanning
CairoMakie.rowsize!(fgl, 1, Relative(0.12))
CairoMakie.rowsize!(fgl, 2, Relative(0.12))
CairoMakie.colgap!(gd, 1, 1)
CairoMakie.rowgap!(gd, 1, 20)
CairoMakie.colsize!(fgl, 1, Relative(0.35))
first_scaffold_axis = gd.content[1].content

ax_land.ylabel = "a. u."
ax_barcodes.ylabel = "Barcode index"

ax_barcodes.xlabel = "Filtration step"
ax_barcodes.xticklabelrotation = 2pi / 6

CairoMakie.xlims!(ax_barcodes, -0.1, x_max_limit + 0.1)
CairoMakie.xlims!(ax_land, -0.01, x_max_limit + 0.1)
ax_land.xticks = (0:max_val)[1:2:end]
ax_barcodes.xticks = (0:max_val)[1:2:end]

f


## ===-===-===-===-===-===-===-===-===-===-===-==-===-===-===-===-===-===-
script_prefix = "1gb3-2"
@info "$(script_prefix): Saving sample matrices"

## ===-===-===-
ENV["GKSwstype"] = "100"

## ===-===-===-
plots_1gb_dir(args...) = plotsdir("section1", script_prefix * "_sample_matrices", args...)

## ===-===-===-
savename_val = "orderd_matrix_topo_properties_and_scaffolds"
pathargs = ("",)

if do_reverse
    pathargs = (pathargs..., "do_reverse=true")
end

## ===-===-===-
landscape_plot_name = plots_1gb_dir(pathargs..., script_prefix * "_$(savename_val).png")
safesave(landscape_plot_name, f)

landscape_plot_name_pdf = plots_1gb_dir(pathargs..., "pdf", script_prefix * "_$(savename_val).pdf")
safesave(landscape_plot_name_pdf, f)

@info "$(script_prefix): Saved."


## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
do_nothing = "ok"
