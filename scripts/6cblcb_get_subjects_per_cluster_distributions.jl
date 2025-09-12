#=
=#

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
# @savename config

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

# config[:total_matrices] = data_info[data_keys[1]] |> get_samples_limit
total_matrices = data_info[data_keys[1]] |> get_samples_limit

force_6cblcb = parsed_args["force_clusters_dist_permutation"]
# if repeat_cycles
# force_6cblcb = true
# end

if any([occursin("fswp", d) for d in data_keys])
    delete!(config, :WASSERSTEIN_DSITANCE)
end

distributions_dictionary, p = produce_or_load(
    datadir("exp_pro", "section6", "distances_distributions", "dim$(selected_dim)"), # path
    config,
    prefix=jld_file_prefix,
    # sort=false,
    force=force_6cblcb,
) do config

    landscapes_distance = populate_dict!(Dict(),
        [clusters_keys],
        # final_structure=
    )
    distances_distribution = populate_dict!(Dict(),
        [clusters_keys],
        final_structure=Float64[]
    )
    # total_matrices = config[:total_matrices]

    # cluster_key = clusters_keys[end]
    cluster_key = "clust95"
    @sync for cluster_key in clusters_keys
        # for cluster_key in clusters_keys
        @info "Cluster key: $(cluster_key)"

        # cluster_landscape_sch = average_landscape_per_cluster[cluster_key]["sch"]
        # cluster_landscape_hc = average_landscape_per_cluster[cluster_key]["hc"]
        cluster_landscape_sch = average_landscape_per_cluster[cluster_key][data_keys[1]]
        cluster_landscape_hc = average_landscape_per_cluster[cluster_key][data_keys[2]]

        # ===-===-===-===-
        # Get the distance- it is computed as the integral of the absolute value of the difference of landscapes
        if isempty(cluster_landscape_hc.land) && isempty(cluster_landscape_sch.land)
            landscapes_distance[cluster_key] = 0
        else
            landscapes_distance[cluster_key] = PersistenceLandscapes.computeDiscanceOfLandscapes(
                cluster_landscape_sch,
                cluster_landscape_hc,
                p_value
            )
        end

        # ===-===-
        total_lands_gr1 = subjects_landscapes_in_clusters[cluster_key][data_keys[1]] |> length
        total_lands_gr2 = subjects_landscapes_in_clusters[cluster_key][data_keys[2]] |> length

        merged_landscapes = vcat(
            [land for (k, land) in subjects_landscapes_in_clusters[cluster_key][data_keys[1]]],
            [land for (k, land) in subjects_landscapes_in_clusters[cluster_key][data_keys[2]]]
        )

        if DO_ZERO_LANDSCAPES
            total_zero_landscapes = 2 * total_matrices - total_lands_gr1 - total_lands_gr2
            merged_landscapes = vcat(
                merged_landscapes,
                [[MyPair(0, 0)] |> PersistenceBarcodes |> PersistenceLandscape for k in 1:total_zero_landscapes]
            )

        end

        total_landscapes = merged_landscapes |> length

        @spawn for shuffle_index in 1:total_shuffles
            # for shuffle_index in 1:total_shuffles
            if shuffle_index % 100 == 0
                @info "\tIteration number= $(shuffle_index)"
                if shuffle_index == total_shuffles
                    @info "\t\tFinishing for cluster = $(cluster_key)"
                end
            end
            shuffling_indices = shuffle(1:total_landscapes)

            if DO_ZERO_LANDSCAPES
                indices_gr1 = shuffling_indices[1:total_matrices]
                indices_gr2 = shuffling_indices[total_matrices+1:end]
            else
                indices_gr1 = shuffling_indices[1:total_lands_gr1]
                indices_gr2 = shuffling_indices[total_lands_gr1+1:end]

                if indices_gr2 |> length != total_lands_gr2
                    ErrorException("Total numbers of cycles left is not same as the second group!") |> throw
                end
            end

            land_vec1 = return_vec_of_landscapes(merged_landscapes, indices_gr1)
            land_vec2 = return_vec_of_landscapes(merged_landscapes, indices_gr2)

            # Get average persistence landscapes for each of the split
            landscape_for_gr1, landscape_for_gr2 = map(x -> x |>
                                                            VectorSpaceOfPersistenceLandscapes |>
                                                            PersistenceLandscapes.real_average,
                [land_vec1, land_vec2]
            )



            if WASSERSTEIN_DSITANCE
                lands_distance =
                    get_wasserstein_form_landscapes(
                        landscape_for_gr1,
                        landscape_for_gr2,
                        p_val=2
                    )
            else
                lands_distance =
                    PersistenceLandscapes.computeDiscanceOfLandscapes(
                        landscape_for_gr1,
                        landscape_for_gr2,
                        p_value
                    )
            end
            push!(distances_distribution[cluster_key], lands_distance)
        end
    end

    Dict(
        "landscapes_distances" => landscapes_distance,
        "distances_distributions" => distances_distribution
    )
end
@info "Keys for the saved distribution df file" keys(distributions_dictionary)

landscapes_distance = distributions_dictionary["landscapes_distances"]
distances_distribution = distributions_dictionary["distances_distributions"]

# Fix clusters keys so that they start from 1 and not 2
if findmin([parse(Int, split(k, "clust")[2]) for (k, v) in distances_distribution])[1] == 2
    new_landscapes_distances = OrderedDict()
    new_distances_distribution = OrderedDict()

    for (k, v) in landscapes_distance
        new_key = "clust$(parse(Int, split(k, "clust")[2])-1)"
        new_landscapes_distances[new_key] = v
    end
    for (k, v) in distances_distribution
        new_key = "clust$(parse(Int, split(k, "clust")[2])-1)"
        new_distances_distribution[new_key] = v
    end

    landscapes_distance = new_landscapes_distances |> sort
    distances_distribution = new_distances_distribution |> sort
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
        # CairoMakie.RGBf(0, 0, 1),
        # CairoMakie.RGBf(142 / 255, 65 / 255, 0 / 255),
        # CairoMakie.RGBf(255 / 255, 165 / 255, 0 / 255),
        # CairoMakie.RGBf(155 / 255, 29 / 255, 137 / 255),
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

critical_values = get_BH_critical_values_form_landscapes(
    landscapes_distance,
    distances_distribution,
)
relevant_cluster_keys = critical_values[critical_values.AcceptedValue, "ClusterKeys"]
selected_clusters = parse.(Int, [split(k, "clust")[2] for k in relevant_cluster_keys])

ordered_permutation =
    [parse(Int, split(c, "clust")[2]) for c in relevant_cluster_keys] |> sortperm
# ===-
# Landscapes and histograms

x_hats = Dict([k => Float64[] for k in data_keys]...)
y_hats = Dict([k => Float64[] for k in data_keys]...)

shifted_start_x = Float64[]
shifted_start_y = Float64[]
shifted_end_x = Float64[]
shifted_end_y = Float64[]
colours_info = Symbol[]

# cl_id = 1
# cluster_key = relevant_cluster_keys[5]
cluster_key = relevant_cluster_keys[ordered_permutation][end-4]

for (cl_id, cluster_key) in relevant_cluster_keys[ordered_permutation] |> enumerate
    cluster_number = parse(Int, split(cluster_key, "clust")[2])

    for (id, selected_key) in data_keys |> enumerate

        pland = average_landscape_per_cluster[cluster_key][selected_key]
        x_vals = [l.first for l in pland.land[1]]
        y_vals = [l.second for l in pland.land[1]]
        x_hat = sum(x_vals .* y_vals) / sum(y_vals)
        y_hat = sum(y_vals .^ 2) / sum(2 * y_vals)

        x_vals = if x_vals == [0.0,]
            @warn "X landscape is empty"
            y_vals
        else
            x_vals
        end
        y_vals = if y_vals == [0.0,]
            @warn "Y landscape is empty"
            x_vals
        else
            y_vals
        end

        push!(x_hats[selected_key], x_hat)
        push!(y_hats[selected_key], y_hat)
    end

    first_key, second_key = if "hc" in data_keys && "sch" in data_keys
        "hc", "sch"
    else
        data_keys[1], data_keys[2]
    end

    arrow_start_x = x_hats[first_key][end]
    arrow_start_y = y_hats[first_key][end]
    arrow_end_x = x_hats[second_key][end]
    arrow_end_y = y_hats[second_key][end]

    push!(shifted_start_x, arrow_start_x .- arrow_start_x)
    push!(shifted_start_y, arrow_start_y .- arrow_start_y)
    push!(shifted_end_x, arrow_end_x .- arrow_start_x)
    push!(shifted_end_y, arrow_end_y .- arrow_start_y)

    # arrow_head_x = x_hats[second_key][end] - x_hats[first_key][end]
    # arrow_head_y = y_hats[second_key][end] - y_hats[first_key][end]
    # arrow_tail_x = x_hats[first_key][end]
    # arrow_tail_y = y_hats[first_key][end]

    # c = get_arrow_colour(arrow_head_x, arrow_head_y)
    # push!(colours_info, c)

end
# Example usage
# x_coords = [1.0, 1.0, 0.0, -1, -1, 0]
# y_coords = [0.0, 1.0, 1.0, 1, 0, -1]
# polar_coords = [cart_to_polar(x, y) for (x, y) in zip(x_coords, y_coords)]
# all_r = [k[1] for k in polar_coords]
# all_theta = ([k[2] for k in polar_coords] ./ π) * 180
# all_theta_rad = [k[2] for k in polar_coords] ./ π
nan_in_x_position = findall(x -> isnan(x), shifted_end_x)
nan_in_y_position = findall(y -> isnan(y), shifted_end_y)

shifted_end_x[nan_in_x_position] .= 0.0
shifted_end_y[nan_in_y_position] .= 0.0


polar_coords = [cart_to_polar(x, y) for (x, y) in zip(shifted_end_x, shifted_end_y)]
all_r = [k[1] for k in polar_coords if !isnan(k[1])]
all_theta = [k[2] for k in polar_coords if !isnan(k[2])]
all_theta_ang = ([k[2] for k in polar_coords if !isnan(k[2])] ./ π) * 180

## ===-===-===-===-===-===-===-===-
# Apply clustering to angles
# test_vectors = [0 0.5 1 0.5 0 -0.5 -1 -2;
#     1 0.5 0 -0.5 -1 -0.5 0 0]
# test_vec_cosine_distance = pairwise(CosineDist(), test_vectors)
# CairoMakie.scatter(test_vectors[1, :], test_vectors[2, :])
scatter_vector_ends = false
if scatter_vector_ends
    max_x = max(shifted_end_x...)
    min_x = min(shifted_end_x...)
    max_y = max(shifted_end_y...)
    min_y = min(shifted_end_y...)

    f = Figure()
    fgl = GridLayout(f[1, 1])
    ax = Axis(fgl[1, 1])
    CairoMakie.xlims!(ax, low=1.1min_x, high=1.1max_x)
    CairoMakie.ylims!(ax, low=1.1min_y, high=1.1max_y)
    for (x, y) in zip(shifted_end_x, shifted_end_y)
        CairoMakie.scatter!(ax, x, y)
        display(f)
    end
end

vector_ends = hcat(
    [x for x in shifted_end_x],
    [y for y in shifted_end_y],)' |> Matrix
vec_cosine_distance = pairwise(CosineDist(), vector_ends)

zeros_in_distances = findall(x -> isnan(x), vec_cosine_distance)
if !isempty(zeros_in_distances)
    @warn "NaN in cosinde distance are replaced with zeros"
end
vec_cosine_distance[zeros_in_distances] .= Inf

link = :average
order = :barjoseph
cos_distance_clustering = hclust(vec_cosine_distance; linkage=link, branchorder=order)

ord = cos_distance_clustering.order
cut_cluster = cutree(cos_distance_clustering, h=0.2)

## ===-
total_landscapes = length(cut_cluster)


sorted_landsacpes_labels =
    [replace(l, "clust" => "Cycle ") for l in
     relevant_cluster_keys[ordered_permutation][ord]
    ]


## ===-===-===-===-===-===-
# Colours setting
arrows_colouring = [
    # PolyElement(color=:grey, strokecolor=:transparent)# , marker='↑'),
    PolyElement(color=:red, strokecolor=:transparent)# , marker='↑'),
    # MarkerElement(color=:black, strokecolor=:transparent)# , marker='↑'),
    # PolyElement(color=:black, strokecolor=:transparent)# , marker='↑'),
    PolyElement(color=:green, strokecolor=:transparent)# , marker='↑'),
]
## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
do_nothing = "ok"
