
## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
using DrWatson
@quickactivate "schtoppaper"

using TopologyPreprocessing
using CairoMakie
using DelimitedFiles
using Statistics

# ===-===-
"MatrixUtils.jl" |> srcdir |> include
"helper_functions.jl" |> srcdir |> include
"PlottingUtils.jl" |> srcdir |> include

# ===-===-

all_cobre_hc = zeros(44, 94, 94)
all_cobre_sch = zeros(44, 94, 94)
some_hcp = zeros(44, 94, 94)

# ===-===-===-===-===-
rawdatadir(args...) = datadir("exp_raw", args...)

# ===-===-===-===-===-
basic_names = ["HC_", "SCH_"]
prefix = "example_"
new_names = "converted_" .* basic_names

# ===-
total_matrices = 44
all_matrices = Dict()
inversed_matrices = Dict()
data_keys = ["hc", "sch"]
for k in data_keys
    all_matrices[k] = Any[]
end

for f in 0:total_matrices-1, k in data_keys
    f_name = prefix * uppercase(string(k)) * "_"
    full_path = rawdatadir("original_weight_matrices", f_name * "$(f)" * ".csv")
    push!(all_matrices[k], readdlm(full_path, ',', Float64, '\n'))
end

for k in data_keys
    inversed_matrices[k] = inverse_elements(all_matrices[k])
end

# ===-===-
max_value_all_matrices = 0
for k in 1:44
    all_cobre_hc[k, :, :] = inversed_matrices["hc"][k]
    all_cobre_sch[k, :, :] = inversed_matrices["sch"][k]

    some_hcp[k, :, :] = readdlm(
        datadir("exp_pro", "original_matrices", "HCP", "original_HCP_$(k).csv"),
        ',',
        Float64,
        '\n'
    )
end
max_value_all_matrices = max(all_cobre_hc[:, :, :]..., all_cobre_sch[:, :, :]..., some_hcp[:, :, :]...)

# ===-===-
total_regions = 94
zeros_counter_hc = []
zeros_counter_sch = []
zeros_counter_hcp = []
for k in 1:44
    local_zeros_counter_hc = 0
    local_zeros_counter_sch = 0
    local_zeros_counter_hcp = 0

    for r in 1:total_regions
        for c in (r+1):total_regions
            if all_cobre_hc[k, r, c] == 0
                local_zeros_counter_hc += 1
            end
            if all_cobre_sch[k, r, c] == 0
                local_zeros_counter_sch += 1
            end
            if some_hcp[k, r, c] == 0
                local_zeros_counter_hcp += 1
            end
        end # c
    end # row
    push!(zeros_counter_hc, local_zeros_counter_hc)
    push!(zeros_counter_sch, local_zeros_counter_sch)
    push!(zeros_counter_hcp, local_zeros_counter_hcp)
end # k


@info "Zeros in HC +/- std: $(mean(zeros_counter_hc)) +/- $(std(zeros_counter_hc))"
@info "Zeros in SCH +/- std: $(mean(zeros_counter_sch)) +/- $(std(zeros_counter_sch))"
@info "Zeros in HCP +/- std: $(mean(zeros_counter_hcp)) +/- $(std(zeros_counter_hcp))"

# all_cobre_hc[coords_above_diagonal...]

# ===-===-

cobre_hc_mat, cobre_sch_mat, hcp_mat = map(
    x -> x |> y -> mean(y, dims=1) |> z -> dropdims(z, dims=1),
    [all_cobre_hc, all_cobre_sch, some_hcp]
)
for k = 1:44
    cobre_hc_mat = all_cobre_hc[k, :, :]
    cobre_sch_mat = all_cobre_sch[k, :, :]
    hcp_mat = some_hcp[k, :, :]

    # ===-===-===-
    # Produce log10 and ordered matrices

    points_per_cm = 28.3465 * 1.6
    plt_width = 17.86 * points_per_cm
    plt_height = 17.86 * points_per_cm
    f = CairoMakie.Figure(size=(plt_width, plt_height,), pt_per_unit=1)
    fgl = GridLayout(f[1, 1])

    CairoMakie.Box(fgl[1, 2], color=:gray90)
    CairoMakie.Label(fgl[1, 2], L"$$ Connectivity matrix", tellwidth=false)
    CairoMakie.Box(fgl[1, 3], color=:gray90)
    CairoMakie.Label(fgl[1, 3], L"$$ Symmetric matrix", tellwidth=false,)
    CairoMakie.Box(fgl[1, 4], color=:gray90)
    CairoMakie.Label(fgl[1, 4], L"$log_{10}$ of symmetric matrix", tellwidth=false,)
    CairoMakie.Box(fgl[1, 5], color=:gray90)
    CairoMakie.Label(fgl[1, 5], L"$$ Ordered matrix", tellwidth=false,)

    CairoMakie.Box(fgl[2, 1], color=:gray90)
    CairoMakie.Label(fgl[2, 1], L"$$ COBRE HC", tellheight=false, rotation=pi / 2,)
    CairoMakie.Box(fgl[3, 1], color=:gray90)
    CairoMakie.Label(fgl[3, 1], L"$$ COBRE SCH", tellheight=false, rotation=pi / 2,)
    CairoMakie.Box(fgl[4, 1], color=:gray90)
    CairoMakie.Label(fgl[4, 1], L"$$ HCP", tellheight=false, rotation=pi / 2,)

    axis_kwargs = (
        aspect=AxisAspect(1),
    )

    mat = cobre_hc_mat
    for (row, mat) in [cobre_hc_mat, cobre_sch_mat, hcp_mat] |> enumerate
        all_zero_connections_indices = findall(x -> x == 0, mat)
        all_zero_connections_indices = vcat(
            [CartesianIndex(n, n) for n in 1:94],
            all_zero_connections_indices
        )

        mat_symmetric = symmetrize_matrix(mat)
        symmetric_mat_with_missing = Matrix{Union{Float64,Missing}}(copy(mat_symmetric))
        all_zero_symmetric_connections_indices = findall(x -> x == 0, mat_symmetric)
        all_zero_symmetric_connections_indices = vcat(
            [CartesianIndex(n, n) for n in 1:94],
            all_zero_symmetric_connections_indices
        )

        symmetric_mat_with_missing[all_zero_symmetric_connections_indices] .= missing


        cobre_hc_mat_with_missing = Matrix{Union{Float64,Missing}}(copy(mat))
        cobre_hc_mat_with_missing[all_zero_connections_indices] .= missing

        cobre_hc_log10 = log10.(mat_symmetric)
        cobre_hc_log10_with_missing = Matrix{Union{Float64,Missing}}(copy(cobre_hc_log10))
        cobre_hc_log10_with_missing[all_zero_symmetric_connections_indices] .= missing

        cobre_hc_ordered = TopologyPreprocessing.get_ordered_matrix(mat_symmetric; assign_same_values=true,
            ordering_start=0
        )
        cobre_hc_ordered_with_missing = Matrix{Union{Float64,Missing}}(copy(cobre_hc_ordered))
        cobre_hc_ordered_with_missing[all_zero_symmetric_connections_indices] .= missing

        ax_cobre_hc = CairoMakie.Axis(fgl[row+1, 1+1][1, 1], title="", aspect=AxisAspect(1))
        ax_cobre_hc_symmetric = CairoMakie.Axis(fgl[row+1, 2+1][1, 1], title="", aspect=AxisAspect(1))
        ax_cobre_hc_log10 = CairoMakie.Axis(fgl[row+1, 3+1][1, 1], title="", aspect=AxisAspect(1))
        ax_cobre_hc_ordered = CairoMakie.Axis(fgl[row+1, 4+1][1, 1], title="", aspect=AxisAspect(1))

        if !issymmetric(cobre_hc_ordered)
            ErrorException("Symmetric matrix is not symmpetric") |> throw
        elseif !issymmetric(cobre_hc_log10)
            ErrorException("Log10 matrix is not symmpetric") |> throw
        elseif !issymmetric(cobre_hc_ordered)
            ErrorException("Ordered matrix is not symmpetric") |> throw
        end

        hm_hc = CairoMakie.heatmap!(
            ax_cobre_hc,
            cobre_hc_mat_with_missing,
            lowclip=:white,
        )
        hm_hc_symmetric = CairoMakie.heatmap!(
            ax_cobre_hc_symmetric,
            symmetric_mat_with_missing,
            lowclip=:white,
        )
        hm_hc_log10 = CairoMakie.heatmap!(
            ax_cobre_hc_log10,
            cobre_hc_log10_with_missing,
            colorrange=(0, log10(max_value_all_matrices)),
            lowclip=:white,
        )
        hm_ord = CairoMakie.heatmap!(
            ax_cobre_hc_ordered,
            cobre_hc_ordered_with_missing,
            lowclip=:white,
        )

        colorbar_kwargs = (
            vertical=false,
            flipaxis=false,
            ticklabelrotation=-pi / 6,
        )
        Colorbar(fgl[row+1, 1+1][2, 1], hm_hc;
            colorbar_kwargs...)
        Colorbar(fgl[row+1, 2+1][2, 1], hm_hc_symmetric;
            colorbar_kwargs...)
        Colorbar(fgl[row+1, 3+1][2, 1], hm_hc_log10;
            colorbar_kwargs...)
        Colorbar(fgl[row+1, 4+1][2, 1], hm_ord;
            tickformat=yticks_formatter,
            colorbar_kwargs...)

        for ax in [ax_cobre_hc, ax_cobre_hc_symmetric, ax_cobre_hc_log10, ax_cobre_hc_ordered]
            hidedecorations!(ax, label=false)
        end
    end # mat

    ## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
    script_prefix = "1k2"
    script_subname = "matrix_comparison_plots"
    @info "$(script_prefix): Saving matrix comparison plots..."

    ## ===-===-===-
    ENV["GKSwstype"] = "100"

    ## ===-===-===-
    plots_1k_dir(args...) = plotsdir("section1", script_prefix * "_$(script_subname)", args...)

    ## ===-===-===-
    joined_keys = "HC-SCH-COBRE"
    f_name = "$(script_subname)_$(joined_keys)_$(k)"

    ## ===-===-===-
    matrix_comparison_fname = plots_1k_dir(script_prefix * "_$(f_name).png")
    safesave(matrix_comparison_fname, f)

    matrix_comparison_fname_pdf = plots_1k_dir("pdf", script_prefix * "_$(f_name).pdf")
    safesave(matrix_comparison_fname_pdf, f)

    @info "$(script_prefix): Saved as" matrix_comparison_fname
end #k

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
do_nothing = "ok"
