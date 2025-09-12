# module COBREConfig

export get_geometric_data_set

## ===-===-===-
import DrWatson: @quickactivate, datadir, srcdir
@quickactivate "schtoppaper"

## ===-===
include(srcdir("BasicConfig.jl"))
include(srcdir("helper_functions.jl"))
include(srcdir("DataSource", "SignalGenerator.jl"))

# using TopologyPreprocessing: generate_random_matrix as get_random_matrix
using TopologyPreprocessing
import DataStructures: OrderedDict
# ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-

# ===-===-===-
function get_geometric_data_set(option::String;
    min_dim::Int=0,
    max_dim::Int=2,
    samples_limiter::Int=1,
    matrix_size::Int=94,
    simulation_len::Int=104
)
    data_info_vector = OrderedDict()
    topo_params = TopoParams(min_dim=min_dim,
        max_dim=max_dim,
        features_list=["norm_barcodes", "bettis"]
    )
    if "1" == option
        # Rand config
        # sym_matrix_source = MatrixGenerator(generate_random_matrix, (total_regions,), 12)

        data_info_vector["rand"] = DataInfo(
            "rand",
            SignalGenerator(
                samples_limiter,
                TopologyPreprocessing.generate_random_matrix,
                (matrix_size,);
                matrix_size=matrix_size,
                total_matrices=simulation_len
            ),
            topo_params;
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
            preprocessing_pipeline=Function[]
        )
    elseif "2" == option
        ## ===-===-
        # Geom config
        cubeR_dim = 3

        data_info_vector["geom_R$(cubeR_dim)"] = DataInfo(
            "geom_R$(cubeR_dim)",
            SignalGenerator(
                samples_limiter,
                get_geometric_matrix,
                (matrix_size, cubeR_dim);
                matrix_size=matrix_size,
                total_matrices=simulation_len
            ),
            topo_params;
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
            preprocessing_pipeline=Function[]
        )

    elseif "3" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [3 5 10 20]

        for R in cubeR_dim
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                SignalGenerator(
                    samples_limiter,
                    get_geometric_matrix,
                    (matrix_size, R);
                    matrix_size=matrix_size,
                    total_matrices=simulation_len
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3a" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [3]

        for R in cubeR_dim
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                SignalGenerator(
                    samples_limiter,
                    get_geometric_matrix,
                    (matrix_size, R);
                    matrix_size=matrix_size,
                    total_matrices=simulation_len
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3b" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [5]

        for R in cubeR_dim
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                SignalGenerator(
                    samples_limiter,
                    get_geometric_matrix,
                    (matrix_size, R);
                    matrix_size=matrix_size,
                    total_matrices=simulation_len
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3c" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [10]

        for R in cubeR_dim
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                SignalGenerator(
                    samples_limiter,
                    get_geometric_matrix,
                    (matrix_size, R);
                    matrix_size=matrix_size,
                    total_matrices=simulation_len
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3d" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [20]

        for R in cubeR_dim
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                SignalGenerator(
                    samples_limiter,
                    get_geometric_matrix,
                    (matrix_size, R);
                    matrix_size=matrix_size,
                    total_matrices=simulation_len
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3ac" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [5 20]

        for R in cubeR_dim
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                SignalGenerator(
                    samples_limiter,
                    get_geometric_matrix,
                    (matrix_size, R);
                    matrix_size=matrix_size,
                    total_matrices=simulation_len
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3e" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [4]

        for R in cubeR_dim
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                SignalGenerator(
                    samples_limiter,
                    get_geometric_matrix,
                    (matrix_size, R);
                    matrix_size=matrix_size,
                    total_matrices=simulation_len
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3f" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [6]

        for R in cubeR_dim
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                SignalGenerator(
                    samples_limiter,
                    get_geometric_matrix,
                    (matrix_size, R);
                    matrix_size=matrix_size,
                    total_matrices=simulation_len
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3g" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [7]

        for R in cubeR_dim
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                SignalGenerator(
                    samples_limiter,
                    get_geometric_matrix,
                    (matrix_size, R);
                    matrix_size=matrix_size,
                    total_matrices=simulation_len
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3h" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [8]

        for R in cubeR_dim
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                SignalGenerator(
                    samples_limiter,
                    get_geometric_matrix,
                    (matrix_size, R);
                    matrix_size=matrix_size,
                    total_matrices=simulation_len
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    else
        @error "Unrecognised option"
    end

    return data_info_vector
end

# end # module
