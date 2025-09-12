#=
=#

import DrWatson: @quickactivate, srcdir, scriptsdir
@quickactivate "schtoppaper"
## ===-===-
"6cb_get_cycles_posets.jl" |> scriptsdir |> include

"0b_convert_regions.jl" |> scriptsdir |> include

"ClusterMerging.jl" |> srcdir |> include
"posets_plotting.jl" |> srcdir |> include
"ArgsParsing.jl" |> srcdir |> include
"helper_functions.jl" |> srcdir |> include
## ===-===-===-===-===-
using PersistenceLandscapes
using Pipe
## ===-===-
# Get 30 clusterings
# 1. what is the clustering that I hve selected
# 2. Get the coloured plot for the selected clustering
# 3. Get the number of the cluster we want to plot.
# 4. Get the average landscape for the selected cluser.

####

parsed_clustering_args = SchiArgPar.parse_clustering_commandline()

@info "6cbl:\tProcessing arguments: " parsed_clustering_args
begin
    @unpack linkage,
    min_clusters,
    cluster_height,
    max_dim,
    min_dim,
    selected_dim,
    total_shuffles,
    p_value,
    zero_landscapes,
    wasserstein_dsitance,
    manual_clustering,
    minimal_height,
    minimal_popularity,
    allow_final_splits,
    max_width,
    min_width,
    do_vector_extension,
    repeat_cycles = parsed_clustering_args
end

MAX_DIM = max_dim
MIN_DIM = resolve_min_dim(min_dim, data_info, data_keys)
linkage = Symbol(linkage)

DO_ZERO_LANDSCAPES = zero_landscapes
WASSERSTEIN_DSITANCE = wasserstein_dsitance

if selected_dim in empty_dimension_warning
    ErrorException("There will be no cycles in selected dimension for given popularity limit. Please change the limit or the selected dimension ") |> throw
end


# ===-===-===-===-===-===-===-===-===-===-===-===-
# Cluster preparation

sdim_index = get_dim_index(MIN_DIM, MAX_DIM, selected_dim)


if do_vector_extension
    total_regions, total_cycles = size(regions_in_cycles_binary_matrix[sdim_index])
    regions_in_cycles = zeros(total_regions + 8, total_cycles)

    # Copy existing info
    regions_in_cycles[1:total_regions, :] = regions_in_cycles_binary_matrix[sdim_index]

    # Add info about yeo networks
    yeo_representation = chceck_yeo_representation(regions_in_cycles_binary_matrix[sdim_index], brain_regions_info)
    regions_in_cycles[(total_regions+1):end, :] = yeo_representation
else
    regions_in_cycles = regions_in_cycles_binary_matrix[sdim_index]
end

total_regions, used_cycles = size(regions_in_cycles)
if repeat_cycles

    popularities = @pipe subject_in_cycle_presence[sdim_index][:, :] |>
                         replace_missing_vals |>
                         binarize_matrix |>
                         sum(_, dims=1) |>
                         vec

    regions_in_cycles = hcat(
        vcat(
            [
                [
                    regions_in_cycles[:, k] for repetition in 1:popularities[k]
                ] for k in 1:used_cycles]
            ...)
        ...)

    subject_presence_in_cycle = subject_in_cycle_presence[sdim_index][:, 1:used_cycles]
    subject_presence_in_cycle = hcat(
        vcat(
            [
                [
                    subject_presence_in_cycle[:, k] for repetition in 1:popularities[k]
                ] for k in 1:used_cycles]
            ...)
        ...)

    ordered_bd_type_split_updated = OrderedDict()
    for local_key in data_keys
        vec1 = ordered_bd_type_split[sdim_index][local_key]
        ordered_bd_type_split_updated[local_key] =
            vcat(
                [
                    [
                        vec1[k] for repetition in 1:popularities[k]
                    ] for k in 1:used_cycles]
                ...)
    end
else
    subject_presence_in_cycle = subject_in_cycle_presence[sdim_index][:, 1:used_cycles]
    ordered_bd_type_split_updated = ordered_bd_type_split[sdim_index]
end


selected_cluster = get_clustering(
    regions_in_cycles;
    selected_linkage=linkage,
    b_order=:r)

if repeat_cycles
    ppl_presence = subject_presence_in_cycle[:, HclustExtended(selected_cluster).order] |> replace_missing_vals

    unique_cycles_info =
        hcat(
            vcat(
                [
                    [
                        ordered_unique_cycles[sdim_index][k] for repetition in 1:popularities[k]
                    ] for k in 1:used_cycles]
                ...)
            ...)
else
    ppl_presence = subject_in_cycle_presence[sdim_index][:, (HclustExtended(selected_cluster)).order] |> replace_missing_vals
    unique_cycles_info = ordered_unique_cycles[sdim_index]
end

if manual_clustering == 0
    clusters_merges = translate_to_merges2(selected_cluster, min_clusters, cluster_height)

elseif manual_clustering > 0
    clusters_merges = split_cluster_n_times(
        HclustExtended(selected_cluster),
        manual_clustering;
        min_height=minimal_height,
        min_popularity=minimal_popularity,
        ppl_presence=ppl_presence,
        allow_final_splits=allow_final_splits,
        max_width=max_width,
        min_width=min_width
    )
else
    ErrorException("manual_clustering have to be 0 or positive, not negative")
end

merged_cluster = translate_merges_to_cutree(clusters_merges, selected_cluster)
merged_cluster .+= 1

unique_clusters_keys = merged_cluster[selected_cluster.order] |> unique
hclust_m = HclustMerged(HclustExtended(selected_cluster), clusters_merges, merged_cluster)

for (i, (k, v)) in hclust_m.clusters |> enumerate
    println("Cluster index $(i), name $(v.name) with range: \t $(v.range)")
end


## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
# Get the persistence data for selected cycles
clusters_keys = ["clust$(c|> get_cluster_split_name)" for (k, c) in hclust_m.clusters]
subject_index = collect(1:(sum([data_info[k] |> get_samples_limit for k in data_keys])))


dictionary_keys =
    [
        clusters_keys,
        data_keys,
    ]
clusters_related_cycles = populate_dict!(Dict(),
    [
        clusters_keys,
    ]
)
barcodes_collection = populate_dict!(Dict(),
    dictionary_keys,
    final_structure=OrderedDict()
)

cluster_cycles_info = populate_dict!(OrderedDict(),
    dictionary_keys,
    final_structure=OrderedDict()
)
subject_related_barcodes = populate_dict!(Dict(),
    [
        clusters_keys,
        data_keys,
        subject_index
    ],
)


for (cname, clust) in hclust_m.clusters
    cluster_index = clust |> get_cluster_split_name
    cluster_key = "clust$(cluster_index)"

    cluster_related_cycles_indices = hclust_m.order[clust.range]

    cluster_cycles = [unique_cycles_info[id] for id in cluster_related_cycles_indices]
    cluster_cycles_len = [length(c[:coords]) for c in cluster_cycles]
    clusters_related_cycles[cluster_key] = cluster_cycles

    ## ===-===-
    # Collect persisntence of cycles

    for (id, cycle_index) in cluster_related_cycles_indices |> enumerate
        for data_key in data_keys
            for k in [barcodes_collection[cluster_key], cluster_cycles_info[cluster_key]]
                k[data_key][cycle_index] = Any[]
            end
        end


        cycle_coords = (cluster_cycles[id][:coords]|>unique)[1]
        for cycle_coords in cluster_cycles[id][:coords] |> unique
            data_key = cycle_coords.key
            pindex = cycle_coords.patient_index
            p_dim = cycle_coords.dim
            cycle_class = cycle_coords.class

            if MIN_DIM == 0
                dim_index = p_dim + 1
            else
                dim_index = p_dim
            end

            # @info patient, local_dim, class
            persistence =
                all_topo_features[data_key][pindex][selected_barcodes][dim_index][cycle_class, :]
            push!(barcodes_collection[cluster_key][data_key][cycle_index],
                persistence
            )

            push!(cluster_cycles_info[cluster_key][data_key][cycle_index],
                ("barcodes" => persistence, "pindex" => pindex, "cycle_class" => cycle_class)
            )

            birth_death =
                all_topo_features[data_key][pindex][selected_barcodes][dim_index][cycle_class, :]

            if !(birth_death in subject_related_barcodes[cluster_key][data_key][pindex])
                push!(subject_related_barcodes[cluster_key][data_key][pindex],
                    birth_death
                )
            end
        end
    end
end
