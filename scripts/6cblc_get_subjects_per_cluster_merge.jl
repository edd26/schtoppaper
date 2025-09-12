#=
=#

import DrWatson: @quickactivate, srcdir, scriptsdir
@quickactivate "schtoppaper"
## ===-===-
"6cbl_extract_barcodes_for_selected_clusters.jl" |> scriptsdir |> include

"ClusterMergePlotting.jl" |> srcdir |> include
## ===-===-===-===-===-
using PersistenceLandscapes

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
# Get the persistence data for selected cycles
@info "6cblc:\tGet the persistence data for selected cycles..."

# Within cluster, get barcodes per subject
# Push to cluster, key, subject

subjects_landscapes_in_clusters = populate_dict!(Dict(), [
        clusters_keys,
        data_keys,
    ],
    final_structure=Dict()
)

for (cname, clust) in hclust_m.clusters
    cluster_index = clust |> get_cluster_split_name
    cluster_key = "clust$(cluster_index)"

    for local_key in data_keys
        local total_matrices = data_info[local_key] |> get_samples_limit

        for pindex in 1:(total_matrices)
            barcodes_vec = subject_related_barcodes[cluster_key][local_key][pindex]
            if !isempty(barcodes_vec)

                barcodes_vec =
                    if barcodes_vec[1][1] == 0.0 && barcodes_vec[1][2] == Inf
                        @warn "Converting infinity barcode to barcode [0,1]"
                        [[0.0, 1.0]]
                    else
                        barcodes_vec
                    end

                subjects_landscapes_in_clusters[cluster_key][local_key][pindex] =
                    @pipe [
                              MyPair(barcode[1], barcode[2])
                              for barcode in barcodes_vec
                          ] |>
                          PersistenceBarcodes(_, sdim_index) |>
                          PersistenceLandscape
            end

        end
    end
end

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
# Create cluster represetation as the average of all subject's lanscapes

average_landscape_per_cluster = populate_dict!(Dict(), [
    clusters_keys,
    data_keys,
],
)
for cluster_key in clusters_keys
    for local_key in data_keys
        vector_space_of_persistence_landscapes = [land for (k, land) in subjects_landscapes_in_clusters[cluster_key][local_key]] |>
                                                 VectorSpaceOfPersistenceLandscapes
        if isempty(vector_space_of_persistence_landscapes.vectOfLand)
            average_landscape = PersistenceLandscape(PersistenceBarcodes([MyPair(0.0, 0.0)]))
        else
            average_landscape = vector_space_of_persistence_landscapes |>
                                PersistenceLandscapes.real_average
        end
        average_landscape_per_cluster[cluster_key][local_key] = average_landscape
    end
end

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
do_nothing = "ok"
