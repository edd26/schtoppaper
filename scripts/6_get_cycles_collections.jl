#=
test Eirene on simple structures (cores and cycloesrep function)

Describe evry of elelemnts with its lifetime

Variables description:
- cycles_collection: raw data from C dictionary for every patient
- unique_cycles: only unique cycles from cycles_collection with coordinates of cycles where it was foudn
- ordered_unique_cycles: the above ordered according to total number of patients for cycle
- nloops_collection_{hc,sch}: unique_cycles with only patients from given group
- filtered_nloops_by_dims_{hc,sch}: above with the cycles in which there are at least X cycles
- filtered_nloops_by_dims: unique_cycles with the cycles in which there are at least X cycles
cycles_type_counters = Vector{Array{Int, 2}}()
ordered_cycles_type_counters

=#

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
import DrWatson: @quickactivate, scriptsdir
@quickactivate "schtoppaper"

using Random
## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
# Including this script exectues all the code included there
"new_init.jl" |> scriptsdir |> include

selected_barcodes = ""
if normalisation_type == to_unique_values::NormalisationType
    selected_barcodes = "norm_barcodes"
elseif normalisation_type == sports::NormalisationType
    selected_barcodes = "sports_norm_barcodes"
elseif normalisation_type == max_elements::NormalisationType
    selected_barcodes = "max_norm_barcodes"
else
    ErrorException("Unrecognised options") |> throw
end
## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-
# using StatsPlots

## ===-===-===-
"SimplicialFaces.jl" |> srcdir |> include
"CyclesUtils.jl" |> srcdir |> include

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
script_prefix = "6"
@info "$(script_prefix):\t Starting cycles analysis..."


# Form the cycles into sets of sets
if !(@isdefined cycles_collection)

    @info "6: Working on cycles_collection..."
    # TODO is the info about the step saved?
    cycles_collection = get_cycles_collection(
        data_info,
        data_keys,
        all_C_data,
        normalisation_type;
        force_production=force_computations,
    )

else
    @info "6:\t Using cycles_collection from workspace."
end

### ===-===-===-===-===-===-===-===-
# Compare the cycles between patietnts and store only unique cycles
# (with informaiton about where it can be found)

# TODO here comes a moment where another config- analysis config- would be useful
if !(@isdefined unique_cycles)
    @info "6:\t Starting unique_cycles extraction..."
    unique_cycles = extract_unique_cycles(cycles_collection; selected_keys=data_keys, min_dim=MIN_DIM, max_dim=MAX_DIM)

else
    @info "6:\t Using unique_cycles from workspace."
end

## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
@info "6:\t Preparing ordered cycles..."
# Subcollections of the unique elements
# Order vectors according to their cardinality
ordered_unique_cycles = copy(unique_cycles);
ordered_unique_cycles_popular = copy(unique_cycles);
ordered_unique_cycles_popular_indices = Dict()# copy(unique_cycles);

limit_popularity = parsed_args["limit_popularity"]
@info "$(script_prefix): Popularity limit: $(limit_popularity)"
empty_dimension_warning = Int[]

pop_threshold = [limit_popularity for k in MIN_DIM:MAX_DIM]
for (dim_index, dim) = MIN_DIM:MAX_DIM |> enumerate
    @info dim
    if isempty(unique_cycles[dim_index])
        ordered_unique_cycles_popular[dim_index] = copy(unique_cycles[dim_index])
        continue
    end
    x = unique_cycles[dim_index][1]

    relevant_keys = data_keys
    ordered_unique_cycles[dim_index] = @pipe [x[:coords] |> unique for x in unique_cycles[dim_index]] |>
                                             [[cc for cc in x if cc.key in relevant_keys] for x in _] |>
                                             [length(x) for x in _] |>
                                             sortperm(_, rev=true) |>
                                             unique_cycles[dim_index][_]

    if findall(x -> x >= pop_threshold[dim_index], [
        x[:coords] |> unique |> length
        for x in ordered_unique_cycles[dim_index]
    ]) |> isempty

        push!(empty_dimension_warning, dim)
        local_threshold = 0
    else
        local_threshold = pop_threshold[dim_index]
    end
    ordered_unique_cycles_popular_indices[dim_index] = popularities =
        findall(x -> x >= local_threshold,
            [
                x[:coords] |> unique |> length
                for x in ordered_unique_cycles[dim_index]
            ])
    ordered_unique_cycles_popular[dim_index] = ordered_unique_cycles[dim_index][popularities]

end

if limit_popularity != 0
    ordered_unique_cycles = ordered_unique_cycles_popular
end

