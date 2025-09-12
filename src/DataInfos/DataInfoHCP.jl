# module COBREConfig
# export data_info

## ===-===-===-
import DrWatson: @quickactivate, datadir
@quickactivate "schtoppaper"

## ===-===
include(srcdir("BasicConfig.jl"))
include(srcdir("helper_functions.jl"))
include(srcdir("DataSource", "MatrixLoader.jl"))


function get_HCP_config(data_label; min_dim::Int=0, max_dim::Int=3, samples_limiter::Int=88)
    features_list = ["barcodes", "bettis", "norm_barcodes", "sports_norm_barcodes", "max_norm_barcodes"]
    ## ===-===-===-
    # Config for a small sample
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices", "HCP"),
        "inversed_HCP_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    small_HCP_config = DataInfo("HCP_small",
        sym_matrix_source,
        topo_params,
        samples_limiter=4,
        ordering_kwargs=(assign_same_values=true, ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-===-
    # Config for 44 subjects
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices", "HCP"),
        "inversed_HCP_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    test_HCP_config = DataInfo("HCP_test",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true, ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-===-
    # Config for 44 subjects, not inversed
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "original_matrices", "HCP"),
        "original_HCP_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    test_HCP_orig_config = DataInfo("HCP_test_orig",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true, ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-===-
    # Config for 44 subjects, not inversed, with detailed matrix ordering
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "original_matrices", "HCP"),
        "original_HCP_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    test_dHCP_orig_config = DataInfo("dHCP_test_orig",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=false, ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )
    ## ===-===-===-
    # All subjects
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices", "HCP"),
        "inversed_HCP_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    main_HCP_config = DataInfo("HCP",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true, ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-===-
    # All subjects
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "original_matrices", "HCP"),
        "original_HCP_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    main_HCP_config_pt1 = DataInfo("HCP_1",
        sym_matrix_source,
        topo_params,
        samples_limiter=samples_limiter,
        ordering_kwargs=(assign_same_values=true, ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    HCP_values_cutoff_config = DataInfo("HCP_1_modified",
        sym_matrix_source,
        topo_params,
        samples_limiter=20,
        ordering_kwargs=(assign_same_values=true, ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, drop_half, reverse_sign],
    )

    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    main_HCP_config_pt1_rev = DataInfo("HCP_1_rev",
        sym_matrix_source,
        topo_params,
        samples_limiter=samples_limiter,
        ordering_kwargs=(assign_same_values=true, ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix],
    )

    sym_matrix_source = MatrixLoader(datadir("exp_pro", "original_matrices", "HCP"),
        "original_HCP_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    main_HCP_config_pt2 = DataInfo("HCP_2",
        sym_matrix_source,
        topo_params,
        samples_limiter=72,
        ordering_kwargs=(assign_same_values=true, ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-
    HCP_data_info = Dict(
        get_identifier(small_HCP_config) => small_HCP_config,
        get_identifier(test_HCP_orig_config) => test_HCP_orig_config,
        get_identifier(test_dHCP_orig_config) => test_dHCP_orig_config,
        get_identifier(test_HCP_config) => test_HCP_config,
        get_identifier(main_HCP_config) => main_HCP_config,
        get_identifier(main_HCP_config_pt1) => main_HCP_config_pt1,
        get_identifier(main_HCP_config_pt2) => main_HCP_config_pt2,
        get_identifier(HCP_values_cutoff_config) => HCP_values_cutoff_config,
        get_identifier(main_HCP_config_pt1_rev) => main_HCP_config_pt1_rev,
    )

    # end # module

    return HCP_data_info[data_label]
end
