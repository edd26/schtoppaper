import DrWatson: @quickactivate, srcdir, scriptsdir
@quickactivate "schtoppaper"
## ===-===-
"6cblc_get_subjects_per_cluster_merge.jl" |> scriptsdir |> include

## ===-===-===-===-===-
using PersistenceLandscapes

using Pipe
import Base.Threads: @threads, @spawn, nthreads, @sync
import CairoMakie: RGBf, PolyElement
import CairoMakie

"DistributionsUtils.jl" |> srcdir |> include
"Centroids.jl" |> srcdir |> include
## ===-===-===-===-===-
@info "6cblcb:\tGet shuffled landscapes's distributions ..."

# ===-===-
get_landscape_from_vector(barcodes_vector; selected_dim=1) = @pipe [MyPair(barcode[1], barcode[2]) for barcode in barcodes_vector] |>
                                                                   PersistenceBarcodes(_, selected_dim) |>
                                                                   PersistenceLandscape

##  Produce or load
jld_file_prefix = "distances_distributions"
joined_keys = join(data_keys, "_")
total_clusters = clusters_keys |> length

if DO_ZERO_LANDSCAPES
    jld_file_prefix, joined_keys = map(
        x -> x * "_zero_padded",
        [jld_file_prefix, joined_keys]
    )
else
    jld_file_prefix, joined_keys = map(
        x -> x * "_not_padded",
        [jld_file_prefix, joined_keys]
    )
end


if manual_clustering > 0
    splits = manual_clustering
    if minimal_height != 0
        config = @dict clusters_keys p_value total_shuffles selected_dim total_clusters joined_keys WASSERSTEIN_DSITANCE splits minimal_height
    else
        config = @dict clusters_keys p_value total_shuffles selected_dim total_clusters joined_keys WASSERSTEIN_DSITANCE splits
    end
else
    config = @dict clusters_keys p_value total_shuffles selected_dim total_clusters joined_keys WASSERSTEIN_DSITANCE
end

if limit_popularity != 0
    config[:minimal_popularity] = limit_popularity
end

if data_set == "humans"
    hcp_key = [k for k in data_keys if occursin("HCP", k)][1]
    config[:total_HCP] = data_info[hcp_key] |> get_samples_limit
end

if ("HCP_test_orig" in data_keys) && do_HCP_trim
    config[:do_HCP_trim] = do_HCP_trim
end

function return_vec_of_landscapes(m_landscapes, indices_gr)
    final_val = m_landscapes[indices_gr]
    if final_val |> isempty
        final_val = [
            PersistenceLandscape([[MyPair(0, 0)]]),
            PersistenceLandscape([[MyPair(0, 0)]])
        ]
    end
    return final_val
end

if repeat_cycles
    config[:repeat_cycles] = repeat_cycles
end

total_matrices = data_info[data_keys[1]] |> get_samples_limit

force_6cblcb = parsed_args["force_clusters_dist_permutation"]

if any([occursin("fswp", d) for d in data_keys])
    delete!(config, :WASSERSTEIN_DSITANCE)
end


data_colours_bank = [
    "hc" => CairoMakie.RGBf(0, 0, 1),
    "sch" => CairoMakie.RGBf(255 / 255, 165 / 255, 0 / 255),
    "HCP_1" => CairoMakie.RGBf(220 / 255, 10 / 255, 20 / 255),
    "COBRE" => CairoMakie.RGBf(120 / 255, 100 / 255, 20 / 255),
] |> OrderedDict

function get_groups_colouring(data_keys, data_colours_bank)
    mri_data_keys = [
        "hc",
        "sch",
        "HCP_1",
        "COBRE",
    ]
    other_colours_bank = [
        cgrad(:BrBG_4, 4, categorical=true)[1:3:4]...,
        cgrad(:grays, max(length(data_keys) - 2, 2), categorical=true)...,
    ]
    others_index = 1

    colors_subjects = []
    for k in data_keys
        if k in mri_data_keys
            push!(colors_subjects, k => data_colours_bank[k])
        else
            push!(colors_subjects, k => other_colours_bank[others_index])
            others_index += 1
        end
    end
    return colors_subjects |> OrderedDict
end

colors_subjects_with_names = get_groups_colouring(data_keys, data_colours_bank)
colors_subjects = [v for (k, v) in colors_subjects_with_names]

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
do_nothing = "ok"
