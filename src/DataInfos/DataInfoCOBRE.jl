# module COBREConfig
# export data_info

## ===-===-===-
import DrWatson: @quickactivate, datadir
@quickactivate "schtoppaper"

## ===-===
include(srcdir("BasicConfig.jl"))
include(srcdir("helper_functions.jl"))
include(srcdir("DataSource", "MatrixLoader.jl"))


function get_COBRE_data_info(data_label; min_dim=0, max_dim=3, samples_limit=44)
    features_list = ["barcodes", "bettis", "norm_barcodes", "sports_norm_barcodes", "max_norm_barcodes"]
    ## ===-===-===-
    # Main hc config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_HC_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    main_hc_config = DataInfo("hc",
        sym_matrix_source,
        topo_params,
        samples_limiter=samples_limit,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-
    # Main sch config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_SCH_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    main_sch_config = DataInfo("sch",
        sym_matrix_source,
        topo_params,
        samples_limiter=samples_limit,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-===-
    # Test hc config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_normalised_HC_",
        "csv"
    )
    topo_params = TopoParams(min_dim=0, max_dim=1, features_list=["norm_barcodes", "bettis"])
    test_hc_config = DataInfo("hc_test",
        sym_matrix_source,
        topo_params,
        samples_limiter=2,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-
    # Test sch config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_normalised_SCH_",
        "csv"
    )
    topo_params = TopoParams(min_dim=0, max_dim=1, features_list=["norm_barcodes", "bettis"])
    test_sch_config = DataInfo("sch_test",
        sym_matrix_source,
        topo_params,
        samples_limiter=5,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-===-
    # Resolution modified hc config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_normalised_HC_",
        "csv"
    )
    topo_params = TopoParams(min_dim=0, max_dim=2, features_list=["norm_barcodes", "bettis"])
    lowres_hc_config2 = DataInfo("hc_lowres2",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    lowres_hc_config4 = DataInfo("hc_lowres4",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    lowres_hc_config10 = DataInfo("hc_lowres10",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )
    lowres_hc_config200 = DataInfo("hc_lowres200",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )
    lowres_hc_config400 = DataInfo("hc_lowres400",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )
    lowres_hc_config800 = DataInfo("hc_lowres800",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )
    ## ===-===-
    # Lowres sch config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_normalised_SCH_",
        "csv"
    )
    topo_params = TopoParams(min_dim=0, max_dim=2, features_list=["norm_barcodes", "bettis"])
    lowres_sch_config2 = DataInfo("sch_lowres2",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    lowres_sch_config4 = DataInfo("sch_lowres4",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )
    lowres_sch_config10 = DataInfo("sch_lowres10",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )
    lowres_sch_config200 = DataInfo("sch_lowres200",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )
    lowres_sch_config400 = DataInfo("sch_lowres400",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )
    lowres_sch_config800 = DataInfo("sch_lowres800",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-===-
    # Detailed hc config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_normalised_HC_",
        "csv"
    )
    topo_params = TopoParams(min_dim=0, max_dim=3, features_list=features_list)
    detailed_hc_config = DataInfo("dhc",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=false,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-
    # Detailed sch config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_normalised_SCH_",
        "csv"
    )
    topo_params = TopoParams(min_dim=0, max_dim=3, features_list=features_list)
    detailed_sch_config = DataInfo("dsch",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(assign_same_values=false,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-===-
    # Not ordered hc config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_normalised_HC_",
        "csv"
    )
    topo_params = TopoParams(min_dim=0, max_dim=3, features_list=features_list)
    nordered_hc_config = DataInfo("hc_not_orded",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(ordering_type=:skip,),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-
    # Not ordered sch config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_normalised_SCH_",
        "csv"
    )
    topo_params = TopoParams(min_dim=0, max_dim=3, features_list=features_list)
    nordered_sch_config = DataInfo("sch_not_orded",
        sym_matrix_source,
        topo_params,
        samples_limiter=44,
        ordering_kwargs=(ordering_type=:skip,),
        preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
    )

    ## ===-===-===-
    # Reversed   hc config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_HC_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    rev_hc_config = DataInfo("rev_hc",
        sym_matrix_source,
        topo_params,
        samples_limiter=samples_limit,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix],
    )

    ## ===-===-
    # Reversed  sch config
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_SCH_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    rev_sch_config = DataInfo("rev_sch",
        sym_matrix_source,
        topo_params,
        samples_limiter=samples_limit,
        ordering_kwargs=(assign_same_values=true,
            ordering_start=0),
        preprocessing_pipeline=Function[symmetrize_matrix,],
    )

    ## ===-===-===-
    # Hc config with zeroed connections replaced with random
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_HC_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    patched_weak_hc_config = DataInfo("patched_weak_hc",
        sym_matrix_source,
        topo_params,
        samples_limiter=samples_limit,
        ordering_kwargs=(assign_same_values=false, ordering_start=0),
        preprocessing_pipeline=Function[randomize_zeros, symmetrize_matrix, reverse_sign,],
    )

    ## ===-===-
    # SCH config with zeroed connections replaced with random
    sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
        "inversed_SCH_",
        "csv"
    )
    topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
    patched_weak_sch_config = DataInfo("patched_weak_sch",
        sym_matrix_source,
        topo_params,
        samples_limiter=samples_limit,
        ordering_kwargs=(assign_same_values=false, ordering_start=0),
        preprocessing_pipeline=Function[randomize_zeros, symmetrize_matrix, reverse_sign,],
    )


    ## ===-===-
    COBRE_data_info = Dict(
        get_identifier(test_hc_config) => test_hc_config,
        get_identifier(test_sch_config) => test_sch_config,
        get_identifier(main_hc_config) => main_hc_config,
        get_identifier(main_sch_config) => main_sch_config,
        get_identifier(lowres_hc_config2) => lowres_hc_config2,
        get_identifier(lowres_sch_config2) => lowres_sch_config2,
        get_identifier(lowres_hc_config4) => lowres_hc_config4,
        get_identifier(lowres_sch_config4) => lowres_sch_config4,
        get_identifier(lowres_hc_config10) => lowres_hc_config10,
        get_identifier(lowres_sch_config10) => lowres_sch_config10,
        get_identifier(lowres_hc_config200) => lowres_hc_config200,
        get_identifier(lowres_sch_config200) => lowres_sch_config200,
        get_identifier(lowres_hc_config400) => lowres_hc_config400,
        get_identifier(lowres_sch_config400) => lowres_sch_config400,
        get_identifier(lowres_hc_config800) => lowres_hc_config800,
        get_identifier(lowres_sch_config800) => lowres_sch_config800,
        get_identifier(detailed_hc_config) => detailed_hc_config,
        get_identifier(detailed_sch_config) => detailed_sch_config,
        get_identifier(nordered_hc_config) => nordered_hc_config,
        get_identifier(nordered_sch_config) => nordered_sch_config,
        get_identifier(rev_hc_config) => rev_hc_config,
        get_identifier(rev_sch_config) => rev_sch_config,
        get_identifier(patched_weak_hc_config) => patched_weak_hc_config,
        get_identifier(patched_weak_sch_config) => patched_weak_sch_config,
    )

    if data_label in keys(COBRE_data_info)
        return COBRE_data_info[data_label]
    else

        sym_matrix_source = MatrixLoader(datadir("exp_pro", "inversed_matrices"),
            data_label,
            "csv"
        )
        topo_params = TopoParams(min_dim=min_dim, max_dim=max_dim, features_list=features_list)
        data_label_config = DataInfo(data_label,
            sym_matrix_source,
            topo_params,
            samples_limiter=44,
            ordering_kwargs=(assign_same_values=true,
                ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )


    end

end
# end # module
