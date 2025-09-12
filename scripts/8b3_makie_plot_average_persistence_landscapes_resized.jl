#=

Plot the average persistence landscapes

=#

import DrWatson: @quickactivate, scriptsdir, @savename
@quickactivate "schtoppaper"
# ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
"8_persistence_landscapes.jl" |> scriptsdir |> include


"ClusterMergePlottingMakie.jl" |> srcdir |> include
# ===-===-===-===-
using CairoMakie

# ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
# Plot the average for both groups
# 1cm <=> 28.3465 pt
# 6.3 * 28.3465

function get_average_pland_plot(dim_range)
    points_per_cm = 28.3465

    plt_width = 17.7 * points_per_cm * 1.5
    plt_height = 10.5 * points_per_cm * 1.5

    f = CairoMakie.Figure(size=(plt_width, plt_height,), pt_per_unit=1)

    fgl = GridLayout(f[1, 1])

    # ===-===- 
    # add lablels 
    row_for_label_shift = 1
    col_for_label_shift = 1

    for (dim_index, local_dim) in dim_range |> enumerate
        Box(fgl[1, dim_index+col_for_label_shift], color=:white)
        Label(fgl[1, dim_index+col_for_label_shift], "Dimension d=$(local_dim)", tellwidth=false)
    end

    for (row_index, label) in [L"Y_{HC}^{d}", L"Y_{SCH}^{d}", L"Y_{HC}^{d}-Y_{SCH}^{d}"] |> enumerate
        Box(fgl[row_index+row_for_label_shift, 1], color=:white)#, color=:gray90)
        Label(fgl[row_index+row_for_label_shift, 1], label, tellheight=false, rotation=pi / 2,)
    end

    for (col, dim_index) in dim_range |> enumerate
        if dim_index == 0
            axis_xlims = (low_x=0 - 0.05, high_x=0.25,)
        else
            axis_xlims = (low_x=0 - 0.1, high_x=1 + 0.1,)
        end

        if dim_index == 0 || dim_index == 1 || dim_index == 2
            axis_ylims = (low_y=0, high_y=0.08)
        elseif dim_index == 3
            axis_ylims = (low_y=0, high_y=0.08)
        end

        for (row, selected_key) in ["hc", "sch"] |> enumerate
            ax_land = CairoMakie.Axis(
                fgl[row+row_for_label_shift, col+col_for_label_shift],
                xtrimspine=true
            )

            ax_land.topspinevisible = false
            ax_land.yticks = 0.00:0.02:0.1
            if dim_index == 0
                ax_land.xticks = 0.0:0.04:1.0
            else
                ax_land.xticks = 0.0:0.2:1.0
            end

            CairoMakie.xlims!(ax_land, axis_xlims...)
            CairoMakie.ylims!(ax_land, axis_ylims...)

            ax_land.rightspinevisible = false
            if col == 1
                ax_land.leftspinevisible = true
                ax_land.ylabel = "a.u."
            else
                hideydecorations!(ax_land, label=true, ticklabels=true,
                    ticks=true, grid=false, minorgrid=false, minorticks=false)
                ax_land.leftspinevisible = false
            end
            pl = average_landscape[selected_key][dim_index]
            plot_persistence_landscape!(ax_land, pl)
        end

        ax_land_diff = CairoMakie.Axis(
            fgl[3+row_for_label_shift, col+col_for_label_shift],
            xlabel="Normalised filtration step",
            xtrimspine=true
        )

        ax_land_diff.yticks = -0.02:0.01:0.02
        if dim_index == 0
            ax_land_diff.xticks = 0.0:0.04:1.0
        else
            ax_land_diff.xticks = 0.0:0.2:1.0
        end

        CairoMakie.xlims!(ax_land_diff, axis_xlims...)
        CairoMakie.ylims!(ax_land_diff, low=-0.02, high=0.02)

        ax_land_diff.rightspinevisible = false
        ax_land_diff.topspinevisible = false
        if col == 1
            ax_land_diff.leftspinevisible = true
            ax_land_diff.ylabel = "a.u."
        else
            hideydecorations!(ax_land_diff, label=true, ticklabels=true,
                ticks=true, grid=false, minorgrid=false, minorticks=false)
            ax_land_diff.leftspinevisible = false
        end

        pl_diff = average_landscape["hc"][dim_index] - average_landscape["sch"][dim_index]
        plot_persistence_landscape!(ax_land_diff, pl_diff)
    end

    return f
end

dim_range = MIN_DIM:MAX_DIM
f_dim0_dim1 = get_average_pland_plot(0:1)
f_dim2_dim3 = get_average_pland_plot(2:3)
# ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
# Save plots
script_prefix = "8b3"

names_to_save = @dict normalisation_type
savename_val = savename(names_to_save, allowedtypes=(DrWatson.default_allowed(names_to_save)..., NormalisationType))

file_prefix = "$(script_prefix)_global_landscapes_with_difference_$(savename_val )"

folder_name = @savename use_outer_layer_only
pathargs = (folder_name,)

safesave(plotsdir("section8", script_prefix * "_global_landscapes_with_difference", pathargs..., file_prefix * "_dim0_dim1" * ".pdf"), f_dim0_dim1)
safesave(plotsdir("section8", script_prefix * "_global_landscapes_with_difference", pathargs..., file_prefix * "_dim2_dim3" * ".pdf"), f_dim2_dim3)
