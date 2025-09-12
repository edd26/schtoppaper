using DrWatson
@quickactivate "schtoppaper"

"1g_landcscapes_samples.jl" |> scriptsdir |> include

using CairoMakie

## ===-===-
rows, _ = size(data_matrices[1])
total_items = ceil(Int, (rows - 1) * rows / 2)
colours_palette = cgrad(:roma, total_items, categorical=true, rev=true);


points_per_cm = 28.3465 * 1.5

plt_width = 19.03 * points_per_cm
plt_height = plt_width * 1.0#(3 / 4)
f = CairoMakie.Figure(size=(plt_width, plt_height,));
fgl = GridLayout(f[1, 1])
matrix_name = ["early-born", "late-born", "long-lived"]

ga = GridLayout(fgl[1, 1])
gb = GridLayout(fgl[2, 1])
gc = GridLayout(fgl[3, 1])

all_grids = [ga, gb, gc]

k = 1
title_sufix = matrix_name[k]
for (k, (fgl_x, title_sufix)) in enumerate(zip(all_grids, matrix_name))

    if do_inverse
        fig_title = "Connectivity matrix,\n $(title_sufix) cycle\n(inv. version)"
    else
        fig_title = "Connectivity matrix,\n $(title_sufix) cycle"
    end

    if do_inverse
        barcode_mat = barcodes[k]
        total_bars, _ = size(barcode_mat)
        land = [[barcode_mat[n, 1], barcode_mat[n, 2]] for n in 1:total_bars]
    else
        land = [refactored_barcodes[k]]
        barcode_mat = barcodes[k]
        total_bars, _ = size(barcode_mat)
        land = [[barcode_mat[n, 1], barcode_mat[n, 2]] for n in 1:total_bars]
    end
    data_mat = data_matrices[k]
    thr = scaffolds_thresholds[k]
    landscape = [[MyPair(l...)] |> PersistenceBarcodes |> PersistenceLandscape for l in land]
    birth_scaffold = data_mat .<= thr.low .&& early_born .!= 0
    death_scaffold = data_mat .< thr.high .&& early_born .!= 0

    alphabet = 'A':'H'
    xticks = yticks = string.(collect(alphabet))

    ax_hmap = Axis(
        fgl_x[1, 1],
        title=fig_title,
        xticks=(1:8, xticks),
        yticks=(1:8, yticks),
        aspect=AxisAspect(1),
    )
    ax_landscape = Axis(
        fgl_x[1, 2],
        title="Persistence landscape,\n $(title_sufix) cycle",
        xtrimspine=true)
    ax_scaffold_birth = Axis(
        fgl_x[1, 3],
        title="Connections at cycle's\nbirth"
    )
    ax_scaffold_death = Axis(
        fgl_x[1, 4],
        title="Connections 1 step\nbefore cycle's death"
    )
    ax_colorbar = fgl_x[1, 5]

    # ===-===-===-
    ax = ax_hmap
    selected_font_size = 10

    # zero elements on diagonal
    hmap_matrix = miss_elements_on_diagonal(data_mat)

    hmap = CairoMakie.heatmap!(ax, hmap_matrix, colormap=colours_palette,)
    for i in 1:8, j in i:8
        if i == j
            continue
        end
        if data_mat[i, j] > 6 && data_mat[i, j] < 23
            txtcolor = :black
        else
            txtcolor = :white
        end
        text!(ax, "$(data_mat[i,j])", position=(i, j), color=txtcolor, align=(:center, :center), fontsize=selected_font_size)
        if i != j
            text!(ax, "$(data_mat[j,i])", position=(j, i), color=txtcolor, align=(:center, :center), fontsize=selected_font_size)
        end
    end
    ax.xticklabelalign = (:right, :center)

    hidexdecorations!(ax, ticks=false, ticklabels=false)
    hideydecorations!(ax, ticks=false, ticklabels=false)
    # ===-====-===-

    # ===-
    for l in landscape
        plot_persistence_landscape!(ax_landscape, l; linewidth=3)
    end

    hidespines!(ax_landscape, :t, :l, :r)

    inf_position = findall(x -> isinf(x), barcodes[1])
    non_inf_barcodes =
        if !isempty(inf_position)
            non_inf_barcodes = copy(barcodes)
            non_inf_barcodes[1] = non_inf_barcodes[1][1:(inf_position[1][1]-1), :]
            non_inf_barcodes
        else
            barcodes
        end

    max_height = ceil(max([max(b...) for b in non_inf_barcodes]...) * √2 / 2)

    CairoMakie.xlims!(ax_landscape, -1, 31)
    CairoMakie.ylims!(ax_landscape, 0, max_height)
    ax_landscape.xlabel = "Filtration step"
    ax_landscape.ylabel = "a.u."

    ax_hmap.yreversed = true
    # ===-

    scaffolds_vec = [birth_scaffold, death_scaffold]
    ax_vector = [ax_scaffold_birth, ax_scaffold_death]

    scaffold_kwargs = (lw=7, alpha=0.8)
    for (scaffold_matrix, ax) in zip(scaffolds_vec, ax_vector)
        CairoMakie.scatter!(ax, x_coords, y_coords)
        append_scaffold(ax, scaffold_matrix, x_coords, y_coords, colours_palette; scaffold_kwargs...)

        hidedecorations!(ax,)
        hidespines!(ax)
    end

    rows, _ = size(data_matrices[1])
    total_items = ceil(Int, (rows - 1) * rows / 2)
    CairoMakie.Colorbar(ax_colorbar,
        limits=(0, total_items),
        colormap=colours_palette,
    )
end

# ===-===-===-===-
# Add plot annotations
for (label, fgl_local) in zip(["a)", "b)", "c)"], [ga, gb, gc])
    Label(fgl_local[1, 1, TopLeft()], label,
        fontsize=18,
        padding=(0, 2, 2, 0),
        halign=:right,
        tellwidth=false
    )
end

CairoMakie.colgap!(fgl, 10)

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
script_prefix = "1g2"
@info "$(script_prefix): Saving sample matrices"

## ===-===-===-
ENV["GKSwstype"] = "100"

## ===-===-===-
plots_1g_dir(args...) = plotsdir("section1", script_prefix * "_sample_matrices", args...)

## ===-===-===-
savename_val = "simple_matrices_landscapes_scaffolds"
pathargs = ("",)

if do_inverse
    pathargs = (pathargs..., "inversed_matrix=true")
end
## ===-===-===-
landscape_plot_name = plots_1g_dir(pathargs..., script_prefix * "_$(savename_val).png")
safesave(landscape_plot_name, f)

landscape_plot_name_pdf = plots_1g_dir(pathargs..., "pdf", script_prefix * "_$(savename_val).pdf")
safesave(landscape_plot_name_pdf, f, pdf_version="1.5")

@info "$(script_prefix): Saved."

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
do_nothing = "ok"
