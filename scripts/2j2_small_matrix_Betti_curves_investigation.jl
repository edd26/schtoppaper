#=


=#

import DrWatson: @quickactivate, scriptsdir, srcdir, @unpack
@quickactivate "schtoppaper"

# ===-===-===-===-===-===-===-===-===-===-===-===-
# "new_init.jl" |> scriptsdir |> include

# ===-===-===-
using CairoMakie
using Base.Iterators: partition
import Base.Threads: @spawn, @sync
using Interpolations
using Pipe
using Eirene

using TopologyPreprocessing
using Random
using DataStructures: OrderedDict

# include(srcdir("helper_functions.jl"))
"ArgsParsingBettiAnalysis.jl" |> srcdir |> include
"OptionStructures.jl" |> srcdir |> include

include(srcdir("DataSource", "MatrixGenerator.jl"))
# include(srcdir("ConfigIndividuals.jl"))
# ===-
# Default configs
datainfodir(args...) = srcdir("DataInfos", args...)
"DataInfoCOBRE.jl" |> datainfodir |> include
"DataInfoHCP.jl" |> datainfodir |> include
# "DataInfoNull.jl" |> srcdir |> include
"DataInfoBOLD.jl" |> datainfodir |> include
# "DataInfoMuldoon2016.jl" |> srcdir |> include
"DataInfoGeometricModel.jl" |> datainfodir |> include
"DataInfoCRT0.jl" |> datainfodir |> include
"DataInfoNull.jl" |> datainfodir |> include

"BettiCurveAnalysis.jl" |> srcdir |> include
# ===-===-===-===-===-===-
script_prefix = "2j2"
script_subname = "Betti_curves_for_reverse_check"
# ===-===-===-
parsed_clustering_args = ArgParseBettis.parse_betti_commandline()
@info "$(script_prefix):\tProcessing arguments: " parsed_clustering_args

begin
    @unpack max_dim,
    min_dim,
    samples_limiter,
    matrix_size,
    regions_side = parsed_clustering_args
end
MAX_DIM = max_dim
MIN_DIM = min_dim

betti_curves_export_path = datadir("exp_pro", "section2", "$(script_prefix)_$(script_subname)")
# ===-===-===-
brain_regions_selection = RegionSelection(regions_side)
total_matrices = samples_limiter

selected_regions =
    if brain_regions_selection == Left
        1:2:2matrix_size
    elseif brain_regions_selection == Right
        2:2:2matrix_size
    elseif brain_regions_selection == Both
        1:matrix_size
    else
        "Region selection not recognised. Possilbe options are $(instances(RegionSelection)), but given: $(brain_regions_selection )" |>
        ErrorException |>
        throw
    end

func_args = Tuple(NaN,)
func_kwargs = NamedTuple()
topo_params = TopoParams(min_dim=MIN_DIM,
    max_dim=MAX_DIM,
    features_list=["norm_barcodes", "bettis"]
)
preproc_funcs = Function[]

# ===-===-
rand_matrices_info = DataInfo(
    "rand",
    MatrixGenerator(
        samples_limiter,
        matrix_size,
        TopologyPreprocessing.generate_random_matrix,
        func_args,
        func_kwargs
    ),
    topo_params;
    samples_limiter=samples_limiter,
    ordering_kwargs=(assign_same_values=true, ordering_start=0),
    # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
    preprocessing_pipeline=preproc_funcs
)

R = 20
geom_matrices_info = DataInfo(
    "geom_R$(R)",
    MatrixGenerator(
        samples_limiter,
        matrix_size,
        get_geometric_matrix,
        (R,),
        func_kwargs
    ),
    topo_params;
    samples_limiter=samples_limiter,
    ordering_kwargs=(assign_same_values=true, ordering_start=0),
    preprocessing_pipeline=Function[]
)

hc_data_info = get_COBRE_data_info("hc", min_dim=MIN_DIM, max_dim=MAX_DIM)
sch_data_info = get_COBRE_data_info("sch", min_dim=MIN_DIM, max_dim=MAX_DIM)

HCP_data_info = get_HCP_config("HCP_1", min_dim=MIN_DIM, max_dim=MAX_DIM)

preproc_funcs_COBRE = hc_data_info |> get_preprocessing
# ===-===-
# Get matrices
rand_symmetric_matrices = rand_matrices_info |> get_matrix_source |> get_matrices
geom_symmetric_matrices = geom_matrices_info |> get_matrix_source |> get_matrices

hc_matrices = (hc_data_info|>get_matrix_source|>get_matrices)[1:total_matrices, selected_regions, selected_regions]
sch_matrices = (sch_data_info|>get_matrix_source|>get_matrices)[1:total_matrices, selected_regions, selected_regions]
HCP_matrices = (HCP_data_info|>get_matrix_source|>get_matrices)[1:total_matrices, selected_regions, selected_regions]

hc_symmetric_matrices = copy(hc_matrices)
sch_symmetric_matrices = copy(sch_matrices)
HCP_symmetric_matrices = copy(HCP_matrices)
for k in 1:total_matrices
    hc_symmetric_matrices[k, :, :] = apply_processing(hc_matrices[k, :, :], preproc_funcs_COBRE)
    sch_symmetric_matrices[k, :, :] = apply_processing(sch_matrices[k, :, :], preproc_funcs_COBRE)
    HCP_symmetric_matrices[k, :, :] = apply_processing(HCP_matrices[k, :, :], preproc_funcs_COBRE)
end

@assert issymmetric(hc_symmetric_matrices[1, :, :])
@assert issymmetric(sch_symmetric_matrices[1, :, :])
@assert issymmetric(HCP_symmetric_matrices[1, :, :])

rand_symmetric_matrices_reversed = .-rand_symmetric_matrices .+ max(rand_symmetric_matrices...)
geom_symmetric_matrices_reversed = .-geom_symmetric_matrices .+ max(geom_symmetric_matrices...)
hc_symmetric_matrices_reversed = .-hc_symmetric_matrices .+ max(hc_symmetric_matrices...)
sch_symmetric_matrices_reversed = .-sch_symmetric_matrices .+ max(sch_symmetric_matrices...)
HCP_symmetric_matrices_reversed = .-HCP_symmetric_matrices .+ max(HCP_symmetric_matrices...)

# ===-===-
# Sanity check
if false
    rand_mat_values = [rand_symmetric_matrices[c] for c in
                       CartesianIndices(rand_symmetric_matrices)
                       if c[2] != c[3]
    ]

    rand_mat_rev_values = [
        rand_symmetric_matrices_reversed[c] for c in
        CartesianIndices(rand_symmetric_matrices_reversed)
        if c[2] != c[3]
    ]

    findmin(rand_mat_values)
    findmax(rand_mat_values)
    findmin(rand_mat_rev_values)
    findmax(rand_mat_rev_values)

end

# ===-===-

data_keys_in_nice_formatting = OrderedDict(
    "rand" => "Random",
    "rev_rand" => "Random, rev.",
    "geomR$(R)" => L"Geom. $R^{%$(R)}$",
    "rev_geomR$(R)" => L"Geom. $R^{%$(R)}$, rev.",
    "hc" => "HC",
    "rev_hc" => "HC, rev.",
    "sch" => "SCH",
    "rev_sch" => "SCH, rev.",
    "HCP" => "HCP",
    "rev_HCP" => "HCP, rev.",
)
data_keys = data_keys_in_nice_formatting |> keys |> collect

# ===-===-
# Get ordering
@info "Working on random"

function produce_bettis(container::Dict)
    # rand_matrices_info
    # max_dim = container[:MAX_DIM]
    # ordered_rand_matrices, rand_matrices_C, rand_params =
    #     get_all_required_data_structures(rand_matrices_info, rand_symmetric_matrices; max_dim=max_dim)
    # ordered_rand_rev_matrices, rand_rev_matrices_C, rand_rev_params =
    #     get_all_required_data_structures(rand_matrices_info, rand_symmetric_matrices_reversed, ; max_dim=max_dim)
    # @strdict ordered_rand_matrices rand_params ordered_rand_rev_matrices rand_rev_params

    max_dim = container[:MAX_DIM]
    local_matrix_slice = container[:matrix_slice]
    data_info = container[:local_data_info]

    ordered_matrices, _, topo_params =
        get_all_required_data_structures(data_info, local_matrix_slice; max_dim=max_dim)
    return Dict("ordered_matrices" => ordered_matrices,
        "topo_params" => topo_params,
    )
end

# ===-===-

data_info_vec = [
    rand_matrices_info, rand_matrices_info,
    rand_matrices_info, rand_matrices_info,
    hc_data_info, hc_data_info,
    sch_data_info, sch_data_info,
    HCP_data_info, HCP_data_info
]
matrices_vec = [
    rand_symmetric_matrices, rand_symmetric_matrices_reversed,
    geom_symmetric_matrices, geom_symmetric_matrices_reversed,
    hc_symmetric_matrices, hc_symmetric_matrices_reversed,
    sch_symmetric_matrices, sch_symmetric_matrices_reversed,
    HCP_symmetric_matrices, HCP_symmetric_matrices_reversed,
]

ordered_rand_matrices, rand_params = Vector{Matrix{Int}}(), Vector{Dict{String,Any}}()
ordered_rand_rev_matrices, rand_rev_params = Vector{Matrix{Int}}(), Vector{Dict{String,Any}}()
ordered_geom_matrices, geom_params = Vector{Matrix{Int}}(), Vector{Dict{String,Any}}()
ordered_geom_rev_matrices, geom_rev_params = Vector{Matrix{Int}}(), Vector{Dict{String,Any}}()
ordered_hc_matrices, hc_params = Vector{Matrix{Int}}(), Vector{Dict{String,Any}}()
ordered_hc_rev_matrices, hc_rev_params = Vector{Matrix{Int}}(), Vector{Dict{String,Any}}()
ordered_sch_matrices, sch_params = Vector{Matrix{Int}}(), Vector{Dict{String,Any}}()
ordered_sch_rev_matrices, sch_rev_params = Vector{Matrix{Int}}(), Vector{Dict{String,Any}}()
ordered_HCP_matrices, HCP_params = Vector{Matrix{Int}}(), Vector{Dict{String,Any}}()
ordered_HCP_rev_matrices, HCP_rev_params = Vector{Matrix{Int}}(), Vector{Dict{String,Any}}()

output_ord_mat_vec = [
    ordered_rand_matrices, ordered_rand_rev_matrices,
    ordered_geom_matrices, ordered_geom_rev_matrices,
    ordered_hc_matrices, ordered_hc_rev_matrices,
    ordered_sch_matrices, ordered_sch_rev_matrices,
    ordered_HCP_matrices, ordered_HCP_rev_matrices,
]
output_topo_features_vec = [
    rand_params, rand_rev_params,
    geom_params, geom_rev_params,
    hc_params, hc_rev_params,
    sch_params, sch_rev_params,
    HCP_params, HCP_rev_params,]

@sync for (local_key, local_data_info, selected_matrices, output_ord_mat, output_topo_features) in zip(data_keys, data_info_vec, matrices_vec, output_ord_mat_vec, output_topo_features_vec)
    @info "Working on $(local_key)"
    @spawn for k in 1:size(selected_matrices, 1)
        # local_key = hcp_key

        matrix_slice = zeros(1, matrix_size, matrix_size)
        matrix_slice[1, :, :] .= selected_matrices[k, :, :]
        container = @dict local_key matrix_size MIN_DIM MAX_DIM matrix_slice local_data_info k

        shuffled_data, s = produce_or_load(
            betti_curves_export_path,
            container, # container
            produce_bettis, # function
            prefix="betti_curves_$(local_key)", # prefix for savename
            tag=false, #github tag
            force=false,
        )

        push!(output_ord_mat, shuffled_data["ordered_matrices"][1])
        push!(output_topo_features, shuffled_data["topo_params"][1])
    end
end

# @assert size(loaded_data["HCP_params"][1]["bettis"], 2) == length(MIN_DIM:MAX_DIM)

# all_C = [
#     rand_matrices_C,
#     rand_rev_matrices_C,
#     geom_matrices_C,
#     geom_rev_matrices_C,
#     hc_matrices_C,
#     hc_rev_matrices_C,
#     sch_matrices,
#     sch_matrices_C,
#     HCP_matrices,
#     HCP_matrices_C
# ]
# ===-===-===-
# 

all_topo_features = OrderedDict(
    "rand" => rand_params,
    "rev_rand" => rand_rev_params,
    "geomR$(R)" => geom_params,
    "rev_geomR$(R)" => geom_rev_params,
    "hc" => hc_params,
    "rev_hc" => hc_rev_params,
    "sch" => sch_params,
    "rev_sch" => sch_rev_params,
    "HCP" => HCP_params,
    "rev_HCP" => HCP_rev_params
)

total_dims = length(MIN_DIM:MAX_DIM)
ylims_per_dim = zeros(total_dims)


k = "rand"
v = all_topo_features[k]
topo_param = v[1]
for (k, v) in all_topo_features
    @info k
    for topo_param in v
        for (dim_index, d) in enumerate(MIN_DIM:MAX_DIM)
            @info "\tDim: $(d)"
            ylims_per_dim[dim_index] = max(
                ylims_per_dim[dim_index],
                max(topo_param["bettis"][:, dim_index]...),
            )
        end
    end
end


# all_topo_features[data_keys[1]][k]["bettis"][:, sdim_index]
# ===-===-===-
betti_colours = TopologyPreprocessing.get_bettis_color_palete(min_dim=MIN_DIM);

total_keys = length(data_keys)
# total_dims = size(all_topo_features[data_keys[1]][1]["bettis"], 2)

# ===-===-===-
base_width = 100
plt_width = base_width * total_keys
plt_height = 110 * total_dims
f = CairoMakie.Figure(size=(plt_width, plt_height));
fgl_main = CairoMakie.GridLayout(f[1, 1])

for (i, d) in enumerate(MIN_DIM:MAX_DIM)
    # CairoMakie.Box(fgl[1, 3], color=:gray90)
    CairoMakie.Label(
        fgl_main[1+i, 1],
        L"D=%$(d)",
        tellheight=false,
        rotation=pi / 2,
        font=:bold
    )
end
for (i, local_key) in enumerate(data_keys)
    # CairoMakie.Box(fgl[1, 3], color=:gray90)
    label = data_keys_in_nice_formatting[local_key]
    CairoMakie.Label(
        fgl_main[1, 1+i],
        L"%$(label)",
        tellwidth=false,
        font=:bold
    )
end
f

fgl = CairoMakie.GridLayout(fgl_main[2:end, 2:end])


total_dims = length(MIN_DIM:MAX_DIM)
total_steps = (matrix_size - 1) * (matrix_size - 2)

low_y = ylims_per_dim .|> log10 .|> ceil .|> y -> max(y, 0.2) |> y -> .-y

for (keys_col, data_key) in data_keys |> enumerate
    betti_axis_per_dim = [
        CairoMakie.Axis(
            fgl[k, keys_col],
            # title="Dimension $(k-1), $(data_key), size=$(matrix_size)"
        ) for k in 1:total_dims]

    for (sdim_index, selected_dim) in MIN_DIM:MAX_DIM |> enumerate
        betti_axis = betti_axis_per_dim[sdim_index]
        col = betti_colours[sdim_index]

        # final_step = round(Int, get_bettis(all_C[keys_col][1], MAX_DIM, min_dim=MIN_DIM)[1][end, 1])
        for k = 1:(total_matrices)

            betti_curve = all_topo_features[data_key][k]["bettis"][:, sdim_index]
            steps = range(0, stop=1, length=max(2, length(betti_curve)))

            CairoMakie.lines!(betti_axis, steps, betti_curve, color=(col, 0.05),)
        end

        all_bettis =
            [all_topo_features[data_key][k]["bettis"][:, sdim_index] for k in 1:total_matrices]
        interpolated_vectors, avg_vector = interpolate_and_average(all_bettis)
        interpolated_steps = range(0, stop=1, length=max(2, length(avg_vector)))

        # Makie.lines!(betti_axis, interpolated_steps, avg_vector, color=(col, 0.9))
        Makie.stairs!(betti_axis, interpolated_steps, avg_vector, color=(col, 0.9), step=:center,)

        # CairoMakie.xlims!(betti_axis, low=0, high=1.1) #highxmax_per_dim[sdim_index])
        CairoMakie.ylims!(betti_axis, low=low_y[sdim_index], high=max(1, 1.1ylims_per_dim[sdim_index]))
        if selected_dim == MAX_DIM
            betti_axis.xlabel = "Normalised\nfiltration step"
        end
        if keys_col == 1
            betti_axis.ylabel = "#cycle"
            hidespines!(betti_axis, :t, :r)
        end
        if keys_col != 1
            hideydecorations!(betti_axis, grid=false)
            hidespines!(betti_axis, :t, :r, :l)
        end


        if selected_dim != MAX_DIM
            hidexdecorations!(betti_axis, grid=false)
        end
        # if selected_dim == MAX_DIM
        # end

    end
    colgap!(fgl, 20)

    if keys_col == 1
        yspace = maximum(tight_yticklabel_spacing!, betti_axis_per_dim)
        for ax in betti_axis_per_dim
            ax.yticklabelspace = yspace
        end
    end
    f

end # data_key


# colgap!(fgl_main, col_dist)
display(f)

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
@info "$(script_prefix): Saving Betti curves..."

## ===-===-===-
ENV["GKSwstype"] = "100"

## ===-===-===-
plots_2j_dir(args...) = plotsdir("section2", script_prefix * "_$(script_subname)", args...)

## ===-===-===-
regions =
    if brain_regions_selection == Left
        "left"
    elseif brain_regions_selection == Right
        "right"
    elseif brain_regions_selection == Both
        "both"
    else
        "Region selection not recognised. Possilbe options are $(instances(RegionSelection)), but given: $(brain_regions_selection )" |>
        ErrorException |>
        throw
    end
dim_range_str = "$(MIN_DIM)-$(MAX_DIM)"
joined_keys = join(data_keys, "-")
pathargs = (
    "dim_range=$(dim_range_str)",
    joined_keys,
)
f_name = "$(script_subname )_$(joined_keys )_dim_range=$(dim_range_str )_size=$(matrix_size)_samplesize=$(total_matrices)_side=$(regions )"
if any(occursin.("prob_thr", data_keys))
    joined_prob = [split(d, "prob_thr")[1] for d in data_keys] |>
                  unique |>
                  (y -> join(y, "_"))
    pathargs = (pathargs...,
        "joined_prob_thr=$(joined_prob )",
    )
    f_name *= "_prob_thr=$(joined_prob)"
else
end


## ===-===-===-
bcurve_plot_name = plots_2j_dir(pathargs..., script_prefix * "_$(f_name).png")
safesave(bcurve_plot_name, f)

bcurve_plot_name_pdf = plots_2j_dir(pathargs..., "pdf", script_prefix * "_$(f_name).pdf")
safesave(bcurve_plot_name_pdf, f)

@info "$(script_prefix): Saved as" bcurve_plot_name

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
do_nothing = "ok"
