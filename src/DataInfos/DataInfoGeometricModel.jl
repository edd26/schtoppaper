# module COBREConfig

export get_geometric_data_set

## ===-===-===-
import DrWatson: @quickactivate, datadir, srcdir
@quickactivate "schtoppaper"

## ===-===
include(srcdir("BasicConfig.jl"))
include(srcdir("helper_functions.jl"))
include(srcdir("DataSource", "MatrixGenerator.jl"))

# using TopologyPreprocessing: generate_random_matrix as get_random_matrix
using TopologyPreprocessing
import DataStructures: OrderedDict
# ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-

# ===-===-===-
function get_geometric_model_data(
    option::String, ;
    min_dim::Int=0,
    max_dim::Int=2,
    samples_limiter::Int=1,
    matrix_size::Int=94,
    func_args::Tuple=Tuple(NaN,),
    func_kwargs::NamedTuple=NamedTuple()
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
            MatrixGenerator(
                samples_limiter,
                matrix_size,
                TopologyPreprocessing.generate_random_matrix,
                func_args,
                func_kwargs
            ),
            topo_params;
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
            preprocessing_pipeline=Function[]
        )
    elseif "1a" == option
        samples_limiter = 22
        # Rand config
        # sym_matrix_source = MatrixGenerator(generate_random_matrix, (total_regions,), 12)
        # func_args = (matrix_size,)
        data_info_vector["rand_1"] =
            DataInfo(
                "rand_1",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        data_info_vector["rand_2"] =
            DataInfo(
                "rand_2",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
    elseif "1b" == option
        # Rand config
        # sym_matrix_source = MatrixGenerator(generate_random_matrix, (total_regions,), 12)
        # func_args = (matrix_size,)
        data_info_vector["rand_1_b30"] =
            DataInfo(
                "rand_1_b30",
                MatrixGenerator(
                    samples_limiter,
                    30,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        data_info_vector["rand_2_b30"] =
            DataInfo(
                "rand_2_b30",
                MatrixGenerator(
                    samples_limiter,
                    30,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
    elseif "1c" == option
        local_matrix_size = 60
        # Rand config
        # sym_matrix_source = MatrixGenerator(generate_random_matrix, (total_regions,), 12)
        # func_args = (matrix_size,)
        data_info_vector["rand_1_b60"] =
            DataInfo(
                "rand_1_b60",
                MatrixGenerator(
                    samples_limiter,
                    local_matrix_size,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        data_info_vector["rand_2_b60"] =
            DataInfo(
                "rand_2_b60",
                MatrixGenerator(
                    samples_limiter,
                    local_matrix_size,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
    elseif "1d" == option
        local_matrix_size = 80
        # Rand config
        # sym_matrix_source = MatrixGenerator(generate_random_matrix, (total_regions,), 12)
        # func_args = (matrix_size,)
        data_info_vector["rand_1_b80"] =
            DataInfo(
                "rand_1_b80",
                MatrixGenerator(
                    samples_limiter,
                    local_matrix_size,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        data_info_vector["rand_2_b80"] =
            DataInfo(
                "rand_2_b80",
                MatrixGenerator(
                    samples_limiter,
                    local_matrix_size,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
    elseif "1e" == option
        samples_limiter = 44
        # Rand config
        # sym_matrix_source = MatrixGenerator(generate_random_matrix, (total_regions,), 12)
        # func_args = (matrix_size,)
        key1 = "rand_$(matrix_size)_1"
        key2 = "rand_$(matrix_size)_2"
        data_info_vector[key1] =
            DataInfo(
                key1,
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        data_info_vector[key2] =
            DataInfo(
                key2,
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
    elseif "1e_rev" == option
        samples_limiter = 44
        # Rand config
        # sym_matrix_source = MatrixGenerator(generate_random_matrix, (total_regions,), 12)
        # func_args = (matrix_size,)
        key1 = "rand_$(matrix_size)_1_rev"
        key2 = "rand_$(matrix_size)_2_rev"
        data_info_vector[key1] =
            DataInfo(
                key1,
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        data_info_vector[key2] =
            DataInfo(
                key2,
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    TopologyPreprocessing.generate_random_matrix,
                    func_args,
                    func_kwargs
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
        func_args = (cubeR_dim,)

        data_info_vector["geom_R$(cubeR_dim)"] = DataInfo(
            "geom_R$(cubeR_dim)",
            MatrixGenerator(
                samples_limiter,
                matrix_size,
                get_geometric_matrix,
                func_args...;
                func_kwargs...
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
            func_args = (R,)
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
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
            func_args = (R,)
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3aa" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [3]

        for R in cubeR_dim
            func_args = (R,)
            data_info_vector["geom_R$(R)_1"] = DataInfo(
                "geom_R$(R)_1",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
            data_info_vector["geom_R$(R)_2"] = DataInfo(
                "geom_R$(R)_2",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
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
            func_args = (R,)
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3bb" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [5]

        for R in cubeR_dim
            func_args = (R,)
            data_info_vector["geom_R$(R)_1"] = DataInfo(
                "geom_R$(R)_1",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
            data_info_vector["geom_R$(R)_2"] = DataInfo(
                "geom_R$(R)_2",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
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
            func_args = (R,)
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3cc" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [10]

        R = cubeR_dim[1]
        func_args = (R,)
        data_info_vector["geom_R$(R)_1"] = DataInfo(
            "geom_R$(R)_1",
            MatrixGenerator(
                samples_limiter,
                matrix_size,
                get_geometric_matrix,
                func_args,
                func_kwargs
            ),
            topo_params;
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
            preprocessing_pipeline=Function[]
        )
        data_info_vector["geom_R$(R)_2"] = DataInfo(
            "geom_R$(R)_2",
            MatrixGenerator(
                samples_limiter,
                matrix_size,
                get_geometric_matrix,
                func_args,
                func_kwargs
            ),
            topo_params;
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
            preprocessing_pipeline=Function[]
        )
    elseif "3d" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [20]

        R = cubeR_dim[1]
        func_args = (R,)
        data_info_vector["geom_R$(R)"] = DataInfo(
            "geom_R$(R)",
            MatrixGenerator(
                samples_limiter,
                matrix_size,
                get_geometric_matrix,
                func_args,
                func_kwargs
            ),
            topo_params;
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
            preprocessing_pipeline=Function[]
        )
    elseif "3dd" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [20]
        R = cubeR_dim[1]

        func_args = (R,)
        data_info_vector["geom_R$(R)_1"] = DataInfo(
            "geom_R$(R)_1",
            MatrixGenerator(
                samples_limiter,
                matrix_size,
                get_geometric_matrix,
                func_args,
                func_kwargs
            ),
            topo_params;
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
            preprocessing_pipeline=Function[]
        )
        data_info_vector["geom_R$(R)_2"] = DataInfo(
            "geom_R$(R)_2",
            MatrixGenerator(
                samples_limiter,
                matrix_size,
                get_geometric_matrix,
                func_args,
                func_kwargs
            ),
            topo_params;
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
            preprocessing_pipeline=Function[]
        )
    elseif "3ac" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [5 20]

        for R in cubeR_dim
            func_args = (R,)
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
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
            func_args = (R,)
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3ee" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [4]

        for R in cubeR_dim
            func_args = (R,)
            data_info_vector["geom_R$(R)_1"] = DataInfo(
                "geom_R$(R)_1",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
            data_info_vector["geom_R$(R)_2"] = DataInfo(
                "geom_R$(R)_2",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
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
            func_args = (R,)
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
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
            func_args = (R,)
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                    func_args,
                    func_kwargs
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
            func_args = (R,)
            data_info_vector["geom_R$(R)"] = DataInfo(
                "geom_R$(R)",
                MatrixGenerator(
                    samples_limiter,
                    matrix_size,
                    get_geometric_matrix,
                ),
                topo_params;
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
                preprocessing_pipeline=Function[]
            )
        end
    elseif "3gg" == option
        ## ===-===-
        # Geom config
        cubeR_dim = [90]
        R = cubeR_dim[1]

        func_args = (R,)
        data_info_vector["geom_R$(R)_1"] = DataInfo(
            "geom_R$(R)_1",
            MatrixGenerator(
                samples_limiter,
                matrix_size,
                get_geometric_matrix,
                func_args,
                func_kwargs
            ),
            topo_params;
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
            preprocessing_pipeline=Function[]
        )
        data_info_vector["geom_R$(R)_2"] = DataInfo(
            "geom_R$(R)_2",
            MatrixGenerator(
                samples_limiter,
                matrix_size,
                get_geometric_matrix,
                func_args,
                func_kwargs
            ),
            topo_params;
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            # preprocessing_pipeline=Function[symmetrize_matrix, z_transform_matrix, reverse_sign]
            preprocessing_pipeline=Function[]
        )
    else
        @error "Unrecognised option"
    end

    return data_info_vector
end

# end # module
