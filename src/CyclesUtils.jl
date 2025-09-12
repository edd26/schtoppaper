
using DrWatson


using Clustering
using Random
using Eirene
using Plots
using Plots.PlotMeasures
import Eirene.Distances: pairwise, Euclidean

struct CycleCoords
    key::String
    patient_index::Int
    dim::Int
    class::Int
end

function get_key(coords::CycleCoords)
    return coords.key
end
function get_index(coords::CycleCoords)
    return coords.patient_index
end
function get_dim(coords::CycleCoords)
    return coords.dim
end
function get_class(coords::CycleCoords)
    return coords.class
end

function get_cycles_nodes_wrap(params_dict::Dict)
    @info "Working..."
    @info "No limiter!"
    @unpack C, max_dim, min_dim, index = params_dict
    results = get_cycles_nodes(C; max_dim=max_dim, min_dim)

    @info "Done."
    return Dict("cycles" => results)
end

"""

Returns a list of nodes for each cycle. Returns a vector for each dimensions in
'C_vec'.
"""
function get_cycles_nodes(C::Dict; max_dim::Int=1, min_dim::Int=1)

    return [[get_maximal_faces(C, class, dim=dim) |> convert_to_set_of_sets
             for class in 1:get_total_classes(C, dim)]
            for dim in min_dim:get_max_dim(C, max_dim)]
end


function get_total_classes(C::Dict, dim::Int)::Int
    @info "\t\t\t\t Next Class"
    length(C["cyclerep"][dim+2])
end

function get_max_dim(C::Dict, dim_limiter::Int)::Int
    @info "\t\t\t Next dim"
    findmax([length(C["cyclerep"]) - 2, dim_limiter])[1]
end

function get_last_C_index(C, matrix_limiter)::Int
    @info "\t\t Next index"
    min(matrix_limiter, length(C)) |> Int
end


function get_maximal_faces(C, class; dim=1, format="vertex x simplex")
    maximal_faces = Eirene.classrep(C, dim=dim, class=class, format=format)
    if isnothing(maximal_faces)
        return Vector{Int}()
    else
        return maximal_faces
    end
end

function convert_to_set_of_sets(some_vector)
    set_of_faces = Set([Set(some_vector[:, k]) for k in 1:size(some_vector, 2)])
end

"""

Returns a list of nodes for each cycle. Returns a vector for each dimensions in
'C_vec'.
"""
function get_cycles_nodes(C_vec::Vector;
    max_sd::Int=-1,
    selected_keys=keys(C_vec),
    total_matrices::Int=-1
)
    cycles_collection = Dict{String,Vector{Vector{Vector{Set{Set{Int}}}}}}()
    for key in selected_keys
        cycles_collection[key] = Vector{Vector{Vector{Set{Set{Int}}}}}()
    end

    for key in selected_keys
        @info "Now on $key"

        total_matrices = length(C_dictionaries[key])


        for t in 1:total_matrices
            @info "\t and $t"
            C = C_dictionaries[key][t]
            max_sd = findmax([length(C["cyclerep"]), max_sd])[1]

            dim_set = Vector{Vector{Set{Set{Int}}}}()
            for sd = 1:max_sd
                class_sets = Vector{Set{Set{Int}}}()
                total_classes = length(C["cyclerep"][sd])

                for calss = 1:total_classes
                    maximal_faces = classrep(C, dim=sd - 2, class=calss, format="vertex x simplex")
                    if !isnothing(maximal_faces)
                        set_of_faces = Set([Set(maximal_faces[:, k]) for k in 1:size(maximal_faces, 2)])
                        push!(class_sets, set_of_faces)
                    end
                end
                push!(dim_set, class_sets)
            end
            push!(cycles_collection[key], dim_set)
        end
    end
end

function convert_to_set_of_sets(some_vector)
    set_of_faces = Set([Set(some_vector[:, k]) for k in 1:size(some_vector, 2)])
end

function merge_coords_with_set!(dim_set, my_set, coordinates_tuple)

    set_location = findall(x -> x[:simplices] == my_set, dim_set)

    # see if set is in collection
    if !isempty(set_location)
        # push coordinate tuple
        for location in set_location
            push!(dim_set[location[1]][:coords], coordinates_tuple)
        end
    else
        # add set to structure
        push!(dim_set, Dict(:simplices => my_set, :coords => CycleCoords[]))

        # push coordinate tuple
        push!(dim_set[end][:coords], coordinates_tuple)
    end #if
end

function extract_unique_cycles(cycles_collection::Union{Dict,OrderedDict};
    selected_keys::Vector{String}=collect(keys(cycles_collection)),
    min_dim::Int=1,
    max_dim::Int=3
)
    unique_cycles = Vector{Vector{Dict{Symbol,Union{Vector{CycleCoords},Set{Set{Int64}}}}}}()
    dim_range = min_dim:max_dim
    for k in dim_range
        push!(unique_cycles, Any[])
    end

    # Collect all cycles that are within cycles generated for sequence of similarity matrices
    @info selected_keys

    for local_key in selected_keys
        @info "Now on $local_key"

        total_matrices = length(cycles_collection[local_key])

        for t in 1:total_matrices
            @info "\t and $t"

            my_set = cycles_collection[local_key][t]
            # Analyse subject's cycles dimension by dimension
            for (dim_index, local_dim) = dim_range |> enumerate
                dim_sets = unique_cycles[dim_index]
                all_classes = size(my_set[dim_index], 1)

                for class = 1:all_classes
                    coordinates_tuple = CycleCoords(local_key, t, local_dim, class)

                    set_location = findall(x -> x[:simplices] == my_set[dim_index][class], dim_sets)

                    # see if set is in collection
                    if !isempty(set_location)
                        # push coordinate tuple
                        for location in set_location
                            push!(dim_sets[location[1]][:coords], coordinates_tuple)
                        end
                    else
                        # add set to structure
                        push!(dim_sets, Dict(:simplices => my_set[dim_index][class], :coords => CycleCoords[]))

                        push!(dim_sets[end][:coords], coordinates_tuple)
                    end #if
                end #class
                unique_cycles[dim_index] = dim_sets
            end # dim_index
        end # t
    end # key

    total_matrices = length(cycles_collection[selected_keys[1]])

    # Aggregate coords of unique cycles
    @info "Aggregate coords of unique cycles"
    for (dim_index, local_dim) = dim_range |> enumerate
        dim_sets = unique_cycles[dim_index]

        for local_key in selected_keys
            @info "Now aggreagating for $local_key, in dimension=$(local_dim)"

            for t in 1:length(cycles_collection[local_key])
                @info "\t and $t"

                my_set = cycles_collection[local_key][t]
                for class = 1:size(my_set[dim_index], 1)
                    merge_coords_with_set!(dim_sets,
                        my_set[dim_index][class],
                        CycleCoords(local_key, t, local_dim, class)
                    )
                end #class
            end # t
        end # key
        unique_cycles[dim_index] = dim_sets
    end # dim

    return unique_cycles
end

# ===== >>>>>

# Filter out single elements
function filter_by_min_elements(unique_cycles; min_elements=1)
    filtered_dims_collection = Any[]
    for k in 1:size(unique_cycles, 1)
        filtered_dim_set = [x for x in unique_cycles[k] if (length(x[:coords]) == 0 || length(x[:coords]) > min_elements)]
        push!(filtered_dims_collection, filtered_dim_set)
    end
    return filtered_dims_collection
end

# filter by patients
function filter_by_patients(unique_cycles; sel_key=:hc, leave_empty=true)
    filtered_dims = Any[]
    for d = 1:length(unique_cycles)
        filtered_dicts = Any[]
        for el = unique_cycles[d]
            selected_coords = [x for x in el[:coords] if sel_key == x.key]
            if !isempty(selected_coords)
                push!(filtered_dicts, Dict(:simplices => el[:simplices], :coords => selected_coords))
            else
                leave_empty && push!(filtered_dicts, Dict(:simplices => el[:simplices], :coords => []))
            end
        end
        push!(filtered_dims, filtered_dicts)
    end
    return filtered_dims
end

function get_unique_vertices_from_cycle(cycle_collection; unique_elements=false)
    regions_vertices = [[x for x in y] for y in cycle_collection]
    if unique_elements
        ab = regions_vertices[1]
        for k in 1:length(regions_vertices)
            ab = vcat(ab, regions_vertices[k])
        end
        return unique(ab)
    else
        return regions_vertices
    end
end


function get_region_names_from_cycle(cycle_collection, regions_labels; unique_elements=false)
    regions_indices = [[x for x in y] for y in cycle_collection]
    regions_names = [[regions_labels[x] for x in y] for y in regions_indices]
    if unique_elements
        if !isempty(regions_names)
            ab = regions_names[1]
            for k in 1:length(regions_indices)
                ab = vcat(ab, regions_names[k])
            end
            return unique(ab)
        else
            return regions_names
        end
    else
        return regions_names
    end
end

function swap_columns!(matrix, target_index, source_index)
    temp_col = matrix[:, target_index]
    matrix[:, target_index] = matrix[:, source_index]
    matrix[:, source_index] = temp_col
    return matrix
end

function swap_elements!(matrix, target_index, source_index)
    temp_col = matrix[target_index]
    matrix[target_index] = matrix[source_index]
    matrix[source_index] = temp_col
    return matrix
end

# ===-===-===-===-
# Posets utils
function define_result_matrix(method, total_regions, total_cycles)
    if method == "lifetimes"
        results_matrix = Array{Union{Float64,Missing},2}(missing, total_regions, total_cycles)
    elseif method == "cardinality"
        results_matrix = Array{Union{Int,Missing},2}(missing, total_regions, total_cycles)
    elseif method == "hamming"
        results_matrix = zeros(Int, total_regions, total_cycles)
    elseif method == "lab"
        results_matrix = zeros(RGB, total_regions, total_cycles)
    elseif method == "rgb"
        results_matrix = zeros(RGB, total_regions, total_cycles)
    end
    return results_matrix
end

function define_result_matrix(method, my_key, total_regions, total_cycles)
    if method == "lifetimes"
        results_matrix = zeros(Float64, total_regions, total_cycles)
    elseif method == "cardinality"
        results_matrix = zeros(Int, total_regions, total_cycles)
    elseif method == "hamming"
        results_matrix = zeros(Int, total_regions, total_cycles)
    elseif method == "rgb"
        results_matrix = zeros(RGB, total_regions, total_cycles)
    end
    return results_matrix
end

# TODO WARNING regions are used in this function as global variable!!!!
function get_regions_vs_cycles(input_data,
    regions_labels::Vector{String};
    method="lifetimes",
    lifetimes=[],
    total_regions::Int=94,
    max_lifetimes=0,
    max_subjects::Int=0,
    do_region_shuffling=Vector{Vector{String}}()
)
    if method == "lifetimes" && isempty(lifetimes)
        throw(Error("For method=$method, lifetimes must be set. Terminating."))
    end
    if method == "rgb" && (max_lifetimes == 0 || max_subjects == 0)
        throw(Error("For method=$method, max_lifetimes and max_subjects must be set. Terminating."))
    end
    if !isempty(do_region_shuffling) # shuffling described above
        @warn "Using regions shuffling"
    end

    total_cycles = length(input_data)
    results_matrix = define_result_matrix(method, total_regions, total_cycles)
    all_indices = collect(1:total_regions)

    for k = 1:total_cycles
        interesing_collection = input_data[k][:simplices]
        if isempty(regions_labels) # This case is to replace gen version of this func
            regions_to_be_located = get_unique_vertices_from_cycle(interesing_collection; unique_elements=true)
        else
            regions_to_be_located = get_region_names_from_cycle(interesing_collection, regions_labels; unique_elements=true)
        end

        for region = regions_to_be_located
            # ===-
            # To do the robustness validation, we allow shuffling of the relation
            # of the rows in the matrix. Since the topology of an ordered matrix
            # and it shuffled* version is the same, we can use the same topological
            # data and only shuffle the ordred of regions. In code, this is equivalent
            # to finding location of the region and adding random intiger from
            # the range (1, total_regions). Or just randomly shuffiling regions_labels
            # before using it to find location (region index).

            if isempty(regions_labels) # This case is to replace gen version of this func
                location = findall(x -> x == region, all_indices)[1]

            elseif !isempty(do_region_shuffling) # shuffling described above
                shuffled_regions = do_region_shuffling[k]
                location = findall(x -> x == region, shuffled_regions)[1]
            else
                location = findall(x -> x == region, regions_labels)[1]
            end
            # ===-


            if method == "lifetimes"
                results_matrix[location, k] = lifetimes[k]

            elseif method == "cardinality"
                results_matrix[location, k] = length(input_data[k][:coords])

            elseif method == "hamming"
                results_matrix[location, k] = 1

            elseif method == "lab"
                a_value = 128 * (lifetimes[k] / max_lifetimes)
                b_value = 128 * (length(input_data[k][:coords]) / max_subjects)
                l_value = 100
                results_matrix[location, k] = Lab(l_value, a_value, b_value)

            elseif method == "rgb"
                r_value = lifetimes[k] / max_lifetimes
                b_value = length(input_data[k][:coords]) / max_subjects
                g_value = 0.0
                results_matrix[location, k] = RGB(r_value, g_value, b_value)
            end
        end
    end

    return results_matrix
end
function get_ppl_vs_cycles(
    ordered_unique_cycles_in_dim,
    total_cycles;
    total_matrices=88,
    all_keys=["hc", "sch"],
    row_shift=[44, 44]
)
    subject_in_cycle_presence = Array{Union{Int,Missing},2}(missing, total_matrices, total_cycles)

    for k in 1:total_cycles
        for coords_tuple = ordered_unique_cycles_in_dim[k][:coords]

            key_index = findfirst(x -> x == coords_tuple.key, all_keys)[1]
            if key_index == 1
                row_index = coords_tuple.patient_index
            else
                row_index = coords_tuple.patient_index + sum(row_shift[1:(key_index-1)])
            end
            subject_in_cycle_presence[row_index, k] = key_index
        end
    end
    return subject_in_cycle_presence
end

function get_ppl_vs_cycles_COBRE(ordered_unique_cycles, total_cycles;
    total_matrices=88
)
    subject_in_cycle_presence = Array{Union{Int,Missing},2}(missing, total_matrices, total_cycles)

    for k in 1:total_cycles
        for coords_tuple = ordered_unique_cycles[k][:coords]
            p_index = coords_tuple.patient_index
            subject_in_cycle_presence[p_index, k] = 1
        end
    end
    return subject_in_cycle_presence
end

# ===-===-===-===-===-===-===-
function get_regions_vs_cycles(input_data,
    regions_labels::Vector{String},
    my_key::String;
    method="lifetimes",
    lifetimes=[],
    total_regions::Int=94,
    max_lifetimes=0,
    max_subjects::Int=0,
    do_region_shuffling=Vector{Vector{String}}()
)
    if method == "lifetimes" && isempty(lifetimes)
        throw(KeyError("For method=$method, lifetimes must be set. Terminating."))
    end
    if method == "rgb" && (max_lifetimes == 0 || max_subjects == 0)
        throw(KeyError("For method=$method, max_lifetimes and max_subjects must be set. Terminating."))
    end
    if !isempty(do_region_shuffling) # shuffling described above
        @warn "Using regions shuffling"
    end

    total_cycles = length(input_data)
    results_matrix = define_result_matrix(method, my_key, total_regions, total_cycles)
    all_indices = collect(1:total_regions)

    for k = 1:total_cycles
        interesing_collection = input_data[k][:simplices]

        if isempty(regions_labels) # This case is to replace gen version of this func
            regions_to_be_located = get_unique_vertices_from_cycle(interesing_collection; unique_elements=true)
        else
            regions_to_be_located = get_region_names_from_cycle(interesing_collection, regions_labels; unique_elements=true)
        end

        for region = regions_to_be_located
            if isempty(regions_labels) # merge with gen function
                location = findall(x -> x == region, all_indices)[1]
            elseif !isempty(do_region_shuffling) # shuffling described above
                shuffled_regions = do_region_shuffling[k]
                location = findall(x -> x == region, shuffled_regions)[1]
            else
                location = findall(x -> x == region, regions_labels)[1]
            end

            # key specific code here:
            all_location_coords = input_data[k][:coords]
            key_specific_coords = [x for x in all_location_coords if x.key == my_key]

            if method == "lifetimes"
                if length(key_specific_coords) > 0
                    results_matrix[location, k] = lifetimes[k]
                end
            elseif method == "cardinality"
                results_matrix[location, k] = length(key_specific_coords)
            elseif method == "hamming"
                if length(key_specific_coords) > 0
                    results_matrix[location, k] = 1
                end
            elseif method == "rgb"
                if length(key_specific_coords) > 0
                    b_value = length(input_data[k][:coords]) / max_subjects
                    r_value = lifetimes[k] / max_lifetimes
                else
                    b_value = 0
                    r_value = 0
                end
                g_value = 0.0
                results_matrix[location, k] = RGB(r_value, g_value, b_value)
            end
        end
    end

    return results_matrix
end

# ===-===-==-

function get_clustering(matrix; distance_metric=Hamming, selected_linkage=:average, b_order=:r, selected_dims=2)
    dm = pairwise(distance_metric(), matrix, dims=selected_dims)

    # normal ordering
    cycles_clustering = hclust(dm, linkage=selected_linkage, branchorder=b_order)
    return cycles_clustering
end

function skip_empty_columns(input_matrix::Matrix)
    total_rows, total_cols = size(input_matrix)
    input_mat_type = typeof(input_matrix[1])

    modified_matix = ([input_matrix[:, k] for k in 1:total_cols if sum(input_matrix[:, k]) != zeros(input_mat_type, 1)[1]])
    new_total_cols = length(modified_matix)

    output_matrix = zeros(input_mat_type, total_rows, new_total_cols)
    for col in 1:new_total_cols
        output_matrix[:, col] = modified_matix[col]
    end
    return output_matrix
end


function plot_clustering(clustering,
    hmap_matrix1,
    hmap_matrix2,
    ord_lifetimes_by_dim;
    max_yticks=88,
    my_palette1=palette([:white, :blue, :orange], max_yticks),
    my_palette2=my_palette1,
    plt_size=(3508, 2480),
    max_xtick=length(ord_lifetimes_by_dim),
    dim_index=1
)
    plot(
        plot(clustering,
            title="Dimension = $(dim)",
            xticks=false),
        heatmap(hmap_matrix1,
            colorbar=false,
            c=my_palette1,
            yticks=(1:10:max_yticks),
            xticks=false,
        ),
        heatmap(hmap_matrix2,
            colorbar=false,
            c=my_palette2,
            yticks=(0:10:max_yticks),
            xticks=false,
        ),
        boxplot(ord_lifetimes_by_dim,
            legend=false,
            borderwidth=0,
            xlims=(1, max_xtick),
            xticks=0:50:max_xtick,
            ylabel="lifetime",
        ),

        layout=grid(4, 1, heights=[0.3, 0.3, 0.3, 0.1]),
        size=plt_size,
        dpi=300,
        leftmargin=10PlotMeasures.mm,
        bottom_margin=10PlotMeasures.mm,
    )
end

function plot_clustering(clustering::Hclust,
    hmap_matrix1::Array,
    hmap_matrix2::Array,
    ppl_matrix::Array,
    unique_nloops::Array,
    ord_lifetimes_by_dim::Vector;
    max_yticks=88,
    my_palette1=palette([cgrad(:roma, 8, categorical=true,)[4], cgrad(:roma, 8, categorical=true,)[7]], 20),
    my_palette2=palette([cgrad(:roma, 8, categorical=true,)[4], cgrad(:roma, 8, categorical=true,)[1]], 20),
    plt_size=(3508, 2480), # vertical A4 in pixels for 300 dpi
    max_xtick=length(ord_lifetimes_by_dim),
    show_clustering=true,
    show_boxplots=true,
    xtick_step=50,
    distributions_included=false,
    dim=1
)
    empty_plot = plot(ticks=nothing, border=:none)

    if show_clustering
        plt1 = plot(clustering,
            xticks=false,
            title="Dimension = $(dim)",
            ylabel="a. u.",
        )
    else
        @warn "Skipping clustering plot"
        plt1 = plot()
    end

    plt2 = heatmap(hmap_matrix1,
        colorbar=false,
        c=my_palette1,
        yticks=(0:10:max_yticks),
        xticks=false,
        ylabel="Region index (popularity)",
    )
    plt3 = heatmap(hmap_matrix2,
        colorbar=false,
        c=my_palette2,
        yticks=(0:10:max_yticks),
        xticks=false,
        ylabel="Region index (lifetimes)",
    )
    plt4 = plot(
        heatmap(ppl_matrix;
            c=palette([:blue, :orange], 2),
            yticks=(0:10:max_yticks),
            ylabel="Patient index",
            xlims=(1, max_xtick),
            legend=false,
            xticks=false
        ),
        groupedbar(unique_nloops;
            bar_position=:stack,
            dpi=300,
            xticks=false,
            xlims=(1, max_xtick),
            bar_width=1,
            yflip=true,
            lw=0,
            legend=false,
            top_margin=-10PlotMeasures.mm,
            border_style=:origin,
            ylabel="total patients"
        ),
        layout=grid(2, 1, heights=[0.7, 0.3]),
    )
    if show_boxplots
        plt5 = boxplot(ord_lifetimes_by_dim,
            legend=false,
            borderwidth=0,
            xlims=(1, max_xtick),
            xticks=0:xtick_step:max_xtick,
            lw=1.5,
            lc=:green,
            markersize=4, # in pixels
            markerstrokewidth=0, # in pixels
            markercolors=:green,
            markeralpha=0.8,
            ylabel="lifetime",
        )
    else
        @info "No Boxplotting, scattering"
        mean_lifetimes = [mean(x) for x in ord_lifetimes_by_dim]
        plt5 = scatter(1:length(mean_lifetimes), mean_lifetimes,
            legend=false,
            borderwidth=0,
            xlims=(1, max_xtick),
            xticks=0:xtick_step:max_xtick,
            lw=0,
            lc=:green,
            markersize=4, # in pixels
            markerstrokewidth=0, # in pixels
            markercolors=:green,
            alpha=0.5,
            ylabel="lifetime",
        )
    end


    if distributions_included
        hist_kwargs = (legend=false,
            orientation=:horizontal,
            lw=0,
            ytick=false,
            yaxis=false,
            left_margin=-20PlotMeasures.mm,
        )
        plt1b = empty_plot
        plt2b = empty_plot

        all_miss = findall(x -> ismissing(x), hmap_matrix2)
        hmap_matrix2z = copy(hmap_matrix2)
        hmap_matrix2z[all_miss] .= 0
        hmap_matrix2z = Float64.(hmap_matrix2z)
        all_non0 = findall(x -> x != 0, hmap_matrix2z)
        hmap_matrix2z[all_non0] .= 1
        plt3b = bar(sum(hmap_matrix2z, dims=2);
            yticks=(0:10:max_yticks),
            ylims=(0, total_regions + 1),
            c=:orange,
            hist_kwargs...)

        ppl_matrixz = copy(ppl_matrix)
        all_miss = findall(x -> ismissing(x), ppl_matrix)
        ppl_matrixz[all_miss] .= 0
        ppl_matrixz = Float64.(ppl_matrixz)
        all_non0 = findall(x -> x != 0, ppl_matrixz)
        ppl_matrixz[all_non0] .= 1
        plt4b = plot(bar(sum(ppl_matrixz, dims=2);
                yticks=(0:10:max_yticks),
                ylims=(0, max_yticks + 1),
                hist_kwargs...),
            empty_plot;
            layout=grid(2, 1, heights=[0.8, 0.2]),
            c=:blue,
            hist_kwargs...)

        ord_lifetimes_by_dimensionz = [mean(ord_lifetimes_by_dim[x]) for x in 1:max_xtick]
        plt5b = histogram(ord_lifetimes_by_dimensionz;
            bins=range(0.0, stop=maximum(ord_lifetimes_by_dimensionz), length=95),
            c=:green,
            hist_kwargs...)

        l = grid(5, 2, heights=[0.2, 0.1, 0.2, 0.3, 0.2], widths=[0.9, 0.1,
            0.9, 0.1,
            0.9, 0.1,
            0.9, 0.1,
            0.9, 0.1,
        ])
    else
        l = grid(5, 1, heights=[0.2, 0.1, 0.2, 0.3, 0.2])
    end
    final_kwargs = (layout=l,
        size=plt_size,
        dpi=300,
        leftmargin=10PlotMeasures.mm,
        bottom_margin=0PlotMeasures.mm,
        thickness_scaling=1.5
    )


    if distributions_included
        return plot(plt1,
            plt1b,
            plt2,
            plt2b,
            plt3,
            plt3b,
            plt4,
            plt4b,
            plt5,
            plt5b;
            final_kwargs...
        )
    else
        return plot(plt1,
            plt2,
            plt3,
            plt4,
            plt5;
            final_kwargs...
        )
    end
end

function plot_clustering_rgb(clustering,
    hmap_matrix1,
    ord_lifetimes_by_dim;
    max_yticks=88,
    plt_size=(3508, 2480),
    max_xtick=length(ord_lifetimes_by_dim),
    dim=1
)
    plot(
        plot(clustering,
            title="Dimension = $(dim)",
            xticks=false),
        heatmap(hmap_matrix1,
            yticks=(0:10:max_yticks),
            xticks=false,
        ),
        boxplot(ord_lifetimes_by_dim,
            legend=false,
            borderwidth=0,
            xlims=(1, max_xtick),
            xticks=0:50:max_xtick,
            ylabel="lifetime",
        ), layout=grid(3, 1, heights=[0.4, 0.4, 0.2]),
        size=plt_size,
        dpi=300,
        leftmargin=10PlotMeasures.mm,
        bottom_margin=10PlotMeasures.mm,
    )

end


function present_legend(; total_data=88, colorspace="rgb")
    legend_x_values = collect(0:total_data)
    legend_y_values = collect(range(0, 1; step=0.05))
    total_y_vals = length(legend_y_values)
    legend_data = zeros(RGB, total_y_vals, total_data)
    if colorspace == "rgb"
        for row in 1:total_y_vals
            for col in 1:total_data
                legend_data[row, col] = RGB(legend_y_values[row], 0, legend_x_values[col])
            end
        end
    elseif colorspace == "lab"
        for row in 1:total_y_vals
            for col in 1:total_data
                legend_data[row, col] = Lab(100,
                    (128 * legend_y_values[row]),
                    (128 * legend_x_values[col]) / total_data
                )
            end
        end
    end
    heatmap(legend_data,
        size=(400, 300),
        xlabel="more cycles ->",
        ylabel="<- more lifetime",
        title="legend for rgb heatmap",
        yflip=true,
    )
end

function replace_zeros_with_missing(matrix::Matrix{T}) where {T<:Real}
    new_matrix = Array{Union{T,Missing},2}(missing, size(matrix))
    non_zeros_location = findall(x -> x != 0, matrix)
    new_matrix[non_zeros_location] .= matrix[non_zeros_location]
    return new_matrix
end


function binarize_matrix(matrix)
    new_matrix = copy(matrix)
    all_miss = findall(x -> ismissing(x), new_matrix)
    if !isempty(all_miss)
        new_matrix[all_miss] .= 0
    end
    all_non0 = findall(x -> x != 0, new_matrix)
    new_matrix[all_non0] .= 1
    return Int.(new_matrix)
end

##
function replace_missing_vals(matrix; replace_value=0)
    all_non_missing = findall(x -> !ismissing(x), matrix)
    new_matrix = zeros(Int, size(matrix))

    for k in all_non_missing
        new_matrix[k] = matrix[k]
    end

    return new_matrix
end

function binarize_matrix(matrix, all_miss)
    new_matrix = copy(matrix)
    new_matrix[all_miss] .= 0
    all_non0 = findall(x -> x != 0, new_matrix)
    new_matrix[all_non0] .= 1
    return Int.(new_matrix)
end


function prepare_for_histograms(mat1, mat2, mat3; nfactor1=1, nfactor2=1, nfactor3=1)
    new_mat1 = binarize_matrix(mat1)
    cvr_sum = sum(new_mat1, dims=1)'

    new_mat2 = binarize_matrix(mat2)
    cvc_sum = sum(new_mat2, dims=1)'

    new_mat3 = copy(mat3)
    new_mat3 = [mean(x) for x in mat3]

    return cvr_sum / nfactor1, cvc_sum / nfactor2, new_mat3 / nfactor3
end

function get_histograms(mat1, mat2, mat3)

    hist_kwargs = (legend=false,
        lw=0,
    )

    plt1 = histogram(mat1;
        c=:orange,
        title="Total regions per cycles histogram",
        xlabel="total regions",
        ylabel="#",
        normalise=true,
        hist_kwargs...)


    plt2 = histogram(mat2;
        c=:orange,
        title="Total patients per cycles histogram",
        xlabel="total patients",
        ylabel="#",
        hist_kwargs...)

    plt3 = histogram(mat3;
        xlims=(0, 1),
        xticks=(0:0.1:1.0),
        c=:green,
        title="Lifetimes histograms",
        xlabel="Lifetime",
        ylabel="#",
        hist_kwargs...)

    return final_plot = plot(plt1, plt2, plt3;
        layout=grid(3, 1, heights=[0.3, 0.3, 0.3]),
        size=(800, 3 * 300),
        dpi=300,
        left_margin=10PlotMeasures.mm,
        bottom_margin=10PlotMeasures.mm
    )
end

function get_cycles_collection(data_info, data_keys, all_C_data, normalisation_type::NormalisationType; force_production=false)
    cycles_collection = populate_dict!(Dict(), [data_keys])

    @info "Producing results..."
    for (index, local_key) in enumerate(data_keys)
        @info "Now on $local_key"
        save_prefix = savename("cycles_nodes", @dict local_key)
        min_dim, max_dim = get_dims_range(data_info[local_key])


        # From C get all cycles with their nodes
        for subject_index in 1:get_samples_limit(data_info[local_key])
            # Get C
            C = copy(all_C_data[local_key][subject_index])
            config = @dict(C,
                max_dim,
                min_dim,
                subject_index,
            )

            config[:normalisation_type] = Symbol(normalisation_type)

            # Get cycles nodes
            loaded_results, p = produce_or_load(
                datadir("6", "cycles_nodes", local_key), # path
                config, # container
                prefix=save_prefix, # prefix for savename
                tag=false, #github tag
                force=force_production,
            ) do config

                @info "Working..."
                @unpack C, max_dim, min_dim, subject_index, normalisation_type = config
                results = get_cycles_nodes(C; max_dim=max_dim, min_dim=min_dim)

                @info "Done."
                loaded_results = Dict("cycles" => results)
            end
        end
    end

    @info "Loading results into scope..."
    for (index, local_key) in enumerate(data_keys)
        @info "Now on $local_key"

        all_cycles = Any[]
        save_prefix = savename("cycles_nodes", @dict local_key)
        min_dim, max_dim = get_dims_range(data_info[local_key])

        # From C get all cycles with their nodes
        for subject_index in 1:get_samples_limit(data_info[local_key])
            # Get C
            C = all_C_data[local_key][subject_index]
            config = @dict(C,
                max_dim,
                min_dim,
                subject_index,
            )

            config[:normalisation_type] = Symbol(normalisation_type)
            # Get cycles nodes
            loaded_results, p = produce_or_load(
                datadir("6", "cycles_nodes", local_key), # path
                config, # container
                prefix=save_prefix, # prefix for savename
                tag=false, #github tag
                force=false,
            ) do config

                ErrorException("This should have been done in the loop above, aborting") |> throw
            end
            cycles_nodes = loaded_results["cycles"]
            push!(all_cycles, cycles_nodes)
        end

        cycles_collection[local_key] = all_cycles# ["all_cycles"]
    end

    return cycles_collection
end # function

# ===-===-
# For 6cblcc

"""
    extract_cycles_of_subject(ordered_unique_cycles, selected_dim, subject_id, data_key)

Returns a vector of `CycleCoords` that contain information on all cycles that are related
to `subject_id`, within `data_key` and `selected_dim`.

"""
function extract_cycles_of_subject(ordered_unique_cycles, selected_dim, subject_id, data_key)# ::Vector{CycleCoords}
    total_cycles = length(ordered_unique_cycles[selected_dim])

    # subject_related_coords = CycleCoords[]
    subject_related_coords = []
    for cycle_id = 1:total_cycles
        for coord in ordered_unique_cycles[selected_dim][cycle_id][:coords] |> unique
            if coord |> get_index == subject_id && coord |> get_key == data_key
                push!(subject_related_coords, cycle_id => coord)
            end # if
        end # coord
    end # cycle_id
    return subject_related_coords
end

"""
    bind_cycle_with_birth(subject_related_coords::Vector{CycleCoords}, all_topo_features::Dict)

Given a vector of `CycleCoords`, returns a vector of paris, where key is the 
cycle id and the value is the birth of the cycle.
"""
function bind_cycle_with_birth(subject_related_cycles, all_topo_features::Dict)::Vector{Pair{Int,Float64}}
    global_classes_with_births = Pair{Int,Float64}[]


    for (global_index, subject_coords) in subject_related_cycles

        data_key, pindex, p_dim, cycle_class =
            map(x -> subject_coords |> x,
                [
                    get_key,
                    get_index,
                    get_dim,
                    get_class])

        birth, death =
            all_topo_features[data_key][pindex]["norm_barcodes"][p_dim][cycle_class, :]

        push!(
            global_classes_with_births,
            global_index => birth
        )
    end
    return global_classes_with_births
end

# ===-===-
