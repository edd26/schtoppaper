
using Combinatorics
# ===-===-===-===-===-===-
# Get all faces for given class

function get_all_simplex_faces(params_dict::Dict)
    @debug "Using DrWatson wrapper"
    @unpack C_data, t, key, barcodes, dims_range = params_dict

    simplex_faces = Dict(:simplex_faces => get_all_simplex_faces(C_data,
        barcodes,
        dims_range))

    return simplex_faces
end

function get_all_simplex_faces2(params_dict::Dict)
    # 2nd version is added because a parameter was added to get_all_simplex_faces
    @debug "Using DrWatson wrapper"
    @unpack C_data, t, key, barcodes, dims_range, combinations = params_dict

    simplex_faces = Dict(:simplex_faces => get_all_simplex_faces(C_data,
        barcodes,
        dims_range,
        do_combinations=combinations))

    return simplex_faces
end

function get_all_simplex_faces(C::Dict, barcodes, dims_range::UnitRange; do_combinations=true)
    """
    get_all_simplex_faces(C::Dict, barcodes, dims_range::UnitRange)

    Takes the simplicial structure analys created in Eirene, stored in 'C', related barcodes
    and a range of dimensions and returns set of sets describing this structure.

    The subset of set are labels of vertices belonging to all faces found in the
    topological structure. Specifically, those are faces of all cycles found in the
    structure at all steps of filtration.

    dims_range must be within range of dimensions used for creating C.
    """
    @debug "Using function"
    set_for_selected_key = Set{Set{Int}}()

    for selected_dim in dims_range
        @debug "\t\tselected_dim $(selected_dim)"
        set_for_selected_dim = Set{Set{Int}}()

        all_classes = length(C["cyclerep"][selected_dim+1]) #+1 because this dictionary starts with dim0

        for selected_class = 1:all_classes # for every cycle in given dimension
            @debug "\t\t\tselected_class $(selected_class)"
            maximal_faces = classrep(C, dim=selected_dim - 1, class=selected_class, format="vertex x simplex")

            # join with selected dim representation
            if do_combinations && !isnothing(maximal_faces)
                # generate all possible faces as Sets{Int}
                lower_faces = get_all_lower_simplex_faces(maximal_faces)

                # Take all common elements for given class and all previous dimensions
                union!(set_for_selected_dim, lower_faces)
            else
                set_of_faces = Set([Set(maximal_faces[:, k]) for k in 1:size(maximal_faces, 2)])
                union!(set_for_selected_dim, set_of_faces)
            end
        end
        # Take all common elements for all previouslu checked dimensions
        union!(set_for_selected_key, set_for_selected_dim)
    end
    return set_for_selected_key
end


# Faces common for all data are called cores

function get_coords_from_cores(cores)
    # convert sets into vectors
    return cores_coords = [[y for y in x] for x in cores if length(x) >= 2]
end

function name_regions_in_cores(edges, regions; do_sorting::Bool=true)
    cores_vertex_names = [[region_label[edge[1]]; region_label[edge[2]]] for edge in edges]

    for edge = 1:length(cores_vertex_names)
        if lowercase(cores_vertex_names[edge][1]) > lowercase(cores_vertex_names[edge][2])
            @info "doing swap for $(edge)"
            cores_vertex_names[edge] = cores_vertex_names[edge][2:-1:1]
        end
    end

    # sort by first element
    do_sorting && sort!(cores_vertex_names)

    return cores_vertex_names
end

function get_vertex_total_connections(cores_vertex_names, region_label)
    total_regions = length(region_label)
    cores_vertex_total_connections = zeros(total_regions)
    # For every vertex check number of connections
    for edge in cores_vertex_names
        index = findall(x -> x == edge[1], region_label)[1]
        cores_vertex_total_connections[index] += 1
        index = findall(x -> x == edge[2], region_label)[1]
        cores_vertex_total_connections[index] += 1
    end
    return cores_vertex_total_connections
end



function get_data_faces(data_keys,
    all_C,
    topo_features,
    data_prefix::String,
    data_path::String;
    dims_range=1:3,
    do_force=false)

    total_matrices = length(all_C[data_keys[1]])
    faces_clean_data = Dict()

    for key = data_keys # for patients group
        @info "Group : $(key)"
        patient_sets = Set{Set{Int}}[]

        for t = 1:total_matrices
            @info "\tPatient: $(t)"

            C_data = all_C[key][t]
            barcodes = topo_features[key][t][:barcodes]
            all_faces1, s = produce_or_load(data_path, # path
                @dict(C_data, t, key, barcodes, dims_range), # container
                get_all_simplex_faces, # function
                prefix="faces_$(data_prefix)", # prefix for savename
                tag=false, #github tag
                force=do_force,
            )

            push!(patient_sets, all_faces1[:simplex_faces])
        end

        faces_clean_data[key] = patient_sets
    end

    return faces_clean_data
end


function get_data_faces2(data_keys,
    all_C,
    topo_features,
    data_prefix::String,
    data_path::String;
    dims_range=1:3,
    do_force=false,
    combinations=true)

    total_matrices = length(all_C[data_keys[1]])
    faces_clean_data = Dict()

    for key = data_keys # for patients group
        @info "Group : $(key)"
        patient_sets = Set{Set{Int}}[]

        for t = 1:total_matrices
            @info "\tPatient: $(t)"

            C_data = all_C[key][t]
            barcodes = topo_features[key][t][:barcodes]
            all_faces1, s = produce_or_load(data_path, # path
                @dict(C_data, t, key, barcodes, dims_range, combinations), # container
                get_all_simplex_faces2, # function
                prefix="faces_$(data_prefix)", # prefix for savename
                tag=false, #github tag
                force=do_force,
            )

            push!(patient_sets, all_faces1[:simplex_faces])
        end

        faces_clean_data[key] = patient_sets
    end

    return faces_clean_data
end


function get_all_lower_simplex_faces(maximal_faces::Array)
    """

    Return a set with all building blocks of the cycle described with simplices
    described in 'maximal_faces'.

    Building bocks are all simplices creating a cycle with all their faces put
    into a set. For example:
        - if a vertex is creating 2 simplicial complexes, in the final set it is
            represented only once
        - if an edge is a face of 2 different simplicial complexes, in the final
            set it is represented only once

    Sample use:
    >>> C = eirene(symmetric_matrix, 'vr')
    >>> # check what are possible classes by studying barcodes
    >>> maximal_faces = classrep(C, dim=2, class=12, format="vertex x simplex")
    >>> class_sets_final = get_all_lower_simplex_faces(maximal_faces)
    """
    class_sets = Set{Int}[]
    combination_size = size(maximal_faces, 1)

    for col = 1:size(maximal_faces, 2)
        building_elements = maximal_faces[:, col]
        for comb_s = combination_size:-1:1

            all_combinations = [Set(k) for k in combinations(building_elements, comb_s)]
            class_sets = [class_sets; all_combinations]
        end
    end

    class_sets_final = Set(class_sets)

    return class_sets_final
end
