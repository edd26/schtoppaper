#=
1. Get regions vs cycles matrix for all data
2. Do the clustering
=#

using DrWatson
@quickactivate "schtoppaper"

# Including this script exectues all the code included there
"6c_get_cycles_lifetimes.jl" |> scriptsdir |> include

##
"CyclesUtils.jl" |> srcdir |> include
"load_coordinates.jl" |> scriptsdir |> include

##
using Clustering
using Eirene.Distances
using Random
using Plots.PlotMeasures

using DataFrames
using StatsBase: median

## ===-===-===-===-===-===-
script_prefix = "6cb"
section_subname = "cycles_popularity_report"
export_dir(args...) = projectdir("logs", "section6", "$(script_prefix)_$(section_subname)", args...)
## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-

function get_dim_index(MIN_DIM, MAX_DIM, selected_dim)
    for (dim_index, dimension) in MIN_DIM:MAX_DIM |> enumerate
        if dimension == selected_dim
            return dim_index
        end
    end
    ErrorException("Failed to find dimension index") |> throw
end

## ===-===-===-===-===-
@info "6cb:\t Getting regions vs cycles matrix for all data"
# 1. Get regions vs cycles matrix for all data
regions_in_cycles_matrix = Any[]
regions_lifetime_per_cycle = Any[]
regions_in_cycles_binary_matrix = Any[]
subject_in_cycle_presence = Any[]

# Get number-region dictionary
file_name = "aal2_center_coords.txt"
number_to_region = get_regions_from_file(file_name)
region_to_number = get_inversed_regions_from_file(file_name)
total_regions = length(region_to_number)

region_label = [number_to_region[k] for k in 1:total_regions]
## ===-===-===-===-
for (dim_index, local_dim) = enumerate(dim_range)
    @info dim_index local_dim
    @info "regions vs cycles"
    push!(regions_in_cycles_matrix,
        get_regions_vs_cycles(
            ordered_unique_cycles[dim_index],
            region_label,
            method="cardinality")
    )
    @info "regions vs cycles lifetimes"
    if isempty(mean_ord_lifetimes[dim_index]) && isempty(ordered_unique_cycles[dim_index])
        regions_vs_cycles_lifetimes_rsult = Float64[]
    else
        regions_vs_cycles_lifetimes_rsult =
            get_regions_vs_cycles(
                ordered_unique_cycles[dim_index],
                region_label,
                method="lifetimes",
                lifetimes=mean_ord_lifetimes[dim_index])
    end
    push!(regions_lifetime_per_cycle, regions_vs_cycles_lifetimes_rsult)

    @info "regions vs cycles for hamming"
    push!(regions_in_cycles_binary_matrix,
        get_regions_vs_cycles(
            ordered_unique_cycles[dim_index],
            region_label,
            method="hamming")
    )

    @info "ppl vs cycles"
    cycles_counter = length(ordered_unique_cycles[dim_index])

    data_keys_popularity = [data_info[k] |> get_samples_limit for k in data_keys]
    push!(subject_in_cycle_presence,
        get_ppl_vs_cycles(ordered_unique_cycles[dim_index],
            # TODO remove cycles counter as a variable- it should be computed from first arg
            cycles_counter;
            total_matrices=sum(data_keys_popularity),
            all_keys=data_keys,
            row_shift=data_keys_popularity
        )
    )
end

## ===-===-
joined_keys = join(data_keys, "_")
prob_thr = ""
if occursin("prob_thr", data_keys[1])
    prob_thr = split(data_keys[1], "prob_thr")[2]
end

f_name = "$(script_prefix)_$(section_subname)_popularity_limit=$(limit_popularity)_$(joined_keys)_max_dim=$(MAX_DIM).txt"
full_file_name = export_dir(joined_keys, f_name)
ispath(export_dir(joined_keys,)) || mkpath(export_dir(joined_keys,))

open(full_file_name, "w") do file
    write(file, "# Data key $(data_keys[1])\n")
    write(file, "# Swap probability $(prob_thr)\n")
end

for (dim_index, selected_dim) in MIN_DIM:MAX_DIM |> enumerate
    @info "Dim: $(selected_dim)"
    mat = subject_in_cycle_presence[dim_index] |> replace_missing_vals
    all_2 = findall(x -> x == 2, mat)
    mat[all_2] .= 1
    popularity_per_cycle = sum(mat, dims=1)
    if isempty(popularity_per_cycle)
        median_result = nothing
        max_val = nothing
        min_val = nothing
    else
        median_result = median(popularity_per_cycle)
        max_val = max(popularity_per_cycle...)
        min_val = min(popularity_per_cycle...)
    end

    ## ===-===-===-===-===-===-
    dim_cycles = ordered_unique_cycles_popular[dim_index]

    dim_cycles_simplices = [d[:simplices] for d in dim_cycles]

    dim_cycles_len = [length(d[:simplices]) for d in dim_cycles]
    dim_cycles_nodes_per_cycle =
        [vcat([[n for n in s] for s in d]...) |> unique |> length for d in dim_cycles_simplices]

    all_lifetimes = [vcat(ordered_lifetimes[dim_index]...)...]

    ordered_birth_times = [[b[1] for b in cycles] for cycles in ordered_bd[dim_index]]
    all_births = [vcat(ordered_birth_times...)...]

    lines = [
        "Total cycles: $(length(popularity_per_cycle))",
        "Mean popularity per cycle is\t $(mean(popularity_per_cycle))",
        "Median popularity per cycle is\t $(median_result)",
        "",
        "Std popularity per cycle is\t $(std(popularity_per_cycle))",
        "Max popularity per cycle is\t $(max_val)",
        "Min popularity per cycle is\t $(min_val)",
        "Popularity >=44 per cycle is\t $(sum(popularity_per_cycle.>=44))",
        "Popularity >=22 per cycle is\t $(sum(popularity_per_cycle.>=22))",
        "Popularity >2 per cycle is\t $(sum(popularity_per_cycle.>2))",
        "Popularity ==2 per cycle is\t $(sum(popularity_per_cycle.==2))",
        "Popularity ==1 per cycle is\t $(sum(popularity_per_cycle.==1))",
        "Popularity >=22 per cycle divided by total cycles is\t $(sum(popularity_per_cycle.>=22)./length(popularity_per_cycle))",
        "",
        "Cycle size is \t$(mean(dim_cycles_len)) +/- $(std(dim_cycles_len))",
        "Nodes per cycle is \t$(mean(dim_cycles_nodes_per_cycle)) +/- $(std(dim_cycles_nodes_per_cycle))",
        "",
        "Mean cycle lifetime:\t $(mean(all_lifetimes)) +/- $(std(all_lifetimes))",
        "Mean cycle birth:\t $(mean(all_births)) +/- $(std(all_births))",
        "",
    ]
    for l in lines
        @info l
    end
    open(full_file_name, "a") do file
        write(file, "## Dimension $(selected_dim)\n")
        for l in lines
            write(file, l * "\n")
        end
    end
end

## ===-===-
do_nothing = "ok"
