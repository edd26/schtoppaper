#=
1. Get lifetimes of all cycles for loop dimension 1:3
2. Compute the mean lifetime for every cycle.
3. Find all cycles which lifetimes is above given value (e.g. 0.1)
4. Find the indices for ordering by total number of patients per cycle.
5. For cycles with mean lifetime above threshold, plot bars stack of class
    cardinality and boxplot of lifetimes.

=#
##

import DrWatson: @quickactivate, scriptsdir
@quickactivate "schtoppaper"

## Including this script exectues all the code included there
"6_get_cycles_collections.jl" |> scriptsdir |> include

## ===-
# using StatsPlots
using Random
using Eirene.Distances

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
@info "6c:\t Starting lifetime extraction..."
# Lifetimes for cycles for every patient check lifetime of the cycle
ordered_lifetimes = Vector{Vector{Vector{Float64}}}()
ordered_bd = Vector{Vector{Vector{Vector{Float64}}}}()
ordered_bd_type_split = Vector{Dict{String,Vector{Vector{Vector{Float64}}}}}()

@info "6c:\t Creating cycles lifetimes..."
for (dim_index, dim) = enumerate(dim_range)
    @info "dim_index: $(dim_index)"
    collection_of_collections_life = Vector{Vector{Float64}}()
    collection_of_collections_bd = Vector{Vector{Vector{Float64}}}()
    collection_of_collections_bd_ts = Dict{String,Vector{Vector{Vector{Float64}}}}()

    for key in data_keys
        collection_of_collections_bd_ts[key] = Vector{Vector{Vector{Float64}}}()
    end

    sets_collection = ordered_unique_cycles[dim_index]
    for data_dict = sets_collection
        lifetimes_collection = Float64[]
        bd_collection = Vector{Vector{Float64}}()
        bd_collection_ts = Dict{String,Vector{Vector{Float64}}}()

        for key in data_keys
            bd_collection_ts[key] = Vector{Vector{Float64}}()
        end

        for coords_tuple = data_dict[:coords] |> unique
            # @info coords_tuple
            key = coords_tuple.key
            patient = coords_tuple.patient_index
            local_dim = coords_tuple.dim
            class = coords_tuple.class

            if MIN_DIM == 0
                dim_index = local_dim + 1
            else
                dim_index = local_dim
            end

            push!(lifetimes_collection,
                all_topo_features[key][patient][selected_barcodes][dim_index][class, 2] -
                all_topo_features[key][patient][selected_barcodes][dim_index][class, 1]
            )
            push!(bd_collection,
                all_topo_features[key][patient][selected_barcodes][dim_index][class, :]
            )
            push!(bd_collection_ts[key],
                all_topo_features[key][patient][selected_barcodes][dim_index][class, :]
            )
        end
        push!(collection_of_collections_life, lifetimes_collection)
        push!(collection_of_collections_bd, bd_collection)
        for key in data_keys

            push!(collection_of_collections_bd_ts[key], bd_collection_ts[key])
        end
    end
    push!(ordered_lifetimes, collection_of_collections_life)
    push!(ordered_bd, collection_of_collections_bd)
    push!(ordered_bd_type_split, collection_of_collections_bd_ts)
end

## ===-===-
@info "6c:\t Creating cycles mean lifetimes..."
mean_ord_lifetimes = Vector{Vector{Float64}}()

for local_dim = 1:length(ordered_lifetimes)
    push!(mean_ord_lifetimes, [mean(x) for x in ordered_lifetimes[local_dim]])
end

##
mean_lifetimes_key_grouped = Dict()

for key in data_keys
    mean_lifetimes_key_grouped[key] = Any[]

    for (dim_index, local_dim) = enumerate(dim_range)
        total_lifetimes = length(mean_ord_lifetimes[dim_index])
        dim_cycles = ordered_unique_cycles[dim_index] |> unique
        key_involved_cycles = Any[]

        for k in 1:total_lifetimes
            key_relevant_data = findall(x -> x.key == key, dim_cycles[k][:coords] |> unique)
            if !isempty(key_relevant_data)

                key_specific_lifetimes = [mean(x) for x in ordered_lifetimes[dim_index][k][key_relevant_data]]
                push!(key_involved_cycles, key_specific_lifetimes)
            end
        end
        push!(mean_lifetimes_key_grouped[key], key_involved_cycles)
    end
end

# END
## ===-===-
do_nothing = "ok"
