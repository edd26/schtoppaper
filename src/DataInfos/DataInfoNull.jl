# module NULLConfig
# export data_info

## ===-===-===-
import DrWatson: @quickactivate, datadir
@quickactivate "schtoppaper"

## ===-===
include(srcdir("BasicConfig.jl"))
include(srcdir("helper_functions.jl"))



function get_null_model_data(
    option::String, ;
    min_dim::Int=0,
    max_dim::Int=2,
    samples_limiter::Int=44,
    matrix_size::Int=94,
    func_args::Tuple=Tuple(NaN,),
    func_kwargs::NamedTuple=NamedTuple(),
    config::String=""
)
    data_info_vector = OrderedDict()
    topo_params = TopoParams(
        min_dim=min_dim,
        max_dim=max_dim,
        features_list=["norm_barcodes", "bettis"]
    )
    if option == "1"
        ## ===-===-===-
        # Main hc config
        sym_matrix_source = MatrixLoader(datadir("exp_raw", "null_model_matrices", "latmio_und_connected"),
            "latmio_und_connected",
            "csv"
        )
        null_model_config_1 = DataInfo("latmio",
            sym_matrix_source,
            topo_params,
            samples_limiter=samples_limiter,
            preprocessing_pipeline=Function[symmetrize_matrix],
        )
        # TODO add info on succesful loading
        ## ===-===-
    elseif option == "2"
        # Main null_model_und_sign config
        sym_matrix_source = MatrixLoader(datadir("exp_raw", "null_model_matrices", "null_model_und_sign"),
            "null_model_und_sign",
            "csv"
        )
        data_info_vector["null_model"] = DataInfo("null_model",
            sym_matrix_source,
            topo_params,
            samples_limiter=samples_limiter,
            preprocessing_pipeline=Function[symmetrize_matrix],
        )
        ## ===-===-
    elseif option == "3"
        # Main randmio_dir_connected config
        sym_matrix_source = MatrixLoader(datadir("exp_raw", "null_model_matrices", "randmio_dir_connected"),
            "randmio_dir_connected",
            "csv"
        )
        data_info_vector["randmio_dir"] = DataInfo(
            "randmio_dir",
            sym_matrix_source,
            topo_params,
            samples_limiter=samples_limiter,
            preprocessing_pipeline=Function[symmetrize_matrix],
        )
        ## ===-===-
    elseif option == "4"
        # Main 
        MatrixLoader(datadir("exp_raw", "null_model_matrices", "randmio_und_connected"),
            "randmio_und_connected",
            "csv"
        )

        data_info_vector["randmio_und"] = DataInfo(
            "randmio_und",
            sym_matrix_source,
            topo_params,
            samples_limiter=samples_limiter,
            preprocessing_pipeline=Function[symmetrize_matrix],
        )

        ## ===-===-
    elseif option == "5"
        # Connectivity based null model, negatve values chagned to zeros
        data_info_vector["nmodel_zero_hc"] = DataInfo(
            "nmodel_zero_hc",
            MatrixLoader(datadir("sims", "connectivity_null_model_zeroed",),
                "simulation_HC_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            preprocessing_pipeline=Function[symmetrize_matrix],
        )

        data_info_vector["nmodel_zero_sch"] = DataInfo(
            "nmodel_zero_sch",
            MatrixLoader(datadir("sims", "connectivity_null_model_zeroed",),
                "simulation_SCH_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            preprocessing_pipeline=Function[symmetrize_matrix],
        )

        ## ===-===-
    elseif option == "6"
        # Connectivity based null model, negatve values chagned to zeros
        data_info_vector["nmodel_abs_hc"] = DataInfo(
            "nmodel_abs_hc",
            MatrixLoader(datadir("sims", "connectivity_null_model_abs",),
                "simulation_abs_HC_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            preprocessing_pipeline=Function[symmetrize_matrix],
        )

        data_info_vector["nmodel_abs_sch"] = DataInfo(
            "nmodel_abs_sch",
            MatrixLoader(datadir("sims", "connectivity_null_model_abs",),
                "simulation_abs_SCH_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            preprocessing_pipeline=Function[symmetrize_matrix],
        )
    elseif option == "7"
        topo_params = TopoParams(
            min_dim=1,
            max_dim=2,
            features_list=["norm_barcodes", "bettis"]
        )
        # Connectivity based null model, negatve values chagned to zeros
        data_info_vector["nmodel_invgauss_zeroed_hc"] = DataInfo(
            "nmodel_invgauss_zeroed_hc",
            MatrixLoader(
                datadir("sims",
                    "connectivity_inversegaussian_model_zeroed",),
                "simulation_HC_",
                "csv"),
            topo_params,
            samples_limiter=15,
            preprocessing_pipeline=Function[symmetrize_matrix],
        )

        data_info_vector["nmodel_invgauss_zeroed_sch"] = DataInfo(
            "nmodel_invgauss_zeroed_sch",
            MatrixLoader(datadir("sims", "connectivity_inversegaussian_model_zeroed",),
                "simulation_SCH_",
                "csv"),
            topo_params,
            samples_limiter=15,
            preprocessing_pipeline=Function[symmetrize_matrix],
        )
    elseif option == "random_rewire_1"
        topo_params = TopoParams(
            min_dim=min_dim,
            max_dim=max_dim,
            features_list=["barcodes", "norm_barcodes", "bettis"]
        )
        # Connectivity based null model, negatve values chagned to zeros
        data_info_vector["random_rewire_1"] = DataInfo(
            "random_rewire_1",
            MatrixLoader(
                datadir("sims",
                    "random_weight_swap",
                    "rewiring_prob=1e-6_total_swaps=100"
                ),
                "null_model_rewiring_prob=1e-6_total_swaps=100_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
    elseif option == "random_rewire_2"
        topo_params = TopoParams(
            min_dim=min_dim,
            max_dim=max_dim,
            features_list=["barcodes", "norm_barcodes", "bettis"]
        )
        # Connectivity based null model, negatve values chagned to zeros
        data_info_vector["random_rewire_2"] = DataInfo(
            "random_rewire_2",
            MatrixLoader(
                datadir("sims",
                    "random_weight_swap",
                    "rewiring_prob=1e-6_total_swaps=1000"
                ),
                "null_model_rewiring_prob=1e-6_total_swaps=1000_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
    elseif option == "random_weight_1"
        topo_params = TopoParams(
            min_dim=min_dim,
            max_dim=max_dim,
            features_list=["barcodes", "norm_barcodes", "bettis"]
        )
        # Connectivity based null model, negatve values chagned to zeros
        data_info_vector["random_weight_hc"] = DataInfo(
            "random_weight_hc",
            MatrixLoader(
                datadir("sims",
                    "COBRE_rnd_connectivity",
                ),
                "hc_rnd_connectivity_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
        data_info_vector["random_weight_sch"] = DataInfo(
            "random_weight_sch",
            MatrixLoader(
                datadir("sims",
                    "COBRE_rnd_connectivity",
                ),
                "sch_rnd_connectivity_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
    elseif option == "random_weight_p_x"
        topo_params = TopoParams(
            min_dim=min_dim,
            max_dim=max_dim,
            features_list=["barcodes", "norm_barcodes", "bettis"]
        )
        # Connectivity based null model, negatve values chagned to zeros

        key_hc = "hc_rnd_connectivity_prob_thr$(config)"
        key_sch = "sch_rnd_connectivity_prob_thr$(config)"
        data_info_vector[key_hc] = DataInfo(
            key_hc,
            MatrixLoader(
                datadir("sims",
                    "COBRE_rnd_connectivity",
                ),
                "hc_rnd_connectivity_prob_thr=$(config)_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
        data_info_vector[key_sch] = DataInfo(
            key_sch,
            MatrixLoader(
                datadir("sims",
                    "COBRE_rnd_connectivity",
                ),
                "sch_rnd_connectivity_prob_thr=$(config)_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
    elseif option == "random_weight_fswp_p_x"
        topo_params = TopoParams(
            min_dim=min_dim,
            max_dim=max_dim,
            features_list=["barcodes", "norm_barcodes", "bettis"]
        )
        # Connectivity based null model, negatve values chagned to zeros

        key_hc = "hc_rnd_connectivity_fswp_prob_thr$(config)"
        key_sch = "sch_rnd_connectivity_fswp_prob_thr$(config)"
        data_info_vector[key_hc] = DataInfo(
            key_hc,
            MatrixLoader(
                datadir("sims",
                    "COBRE_rnd_connectivity",
                ),
                "hc_rnd_connectivity_fswp_prob_thr=$(config)_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
        data_info_vector[key_sch] = DataInfo(
            key_sch,
            MatrixLoader(
                datadir("sims",
                    "COBRE_rnd_connectivity",
                ),
                "sch_rnd_connectivity_fswp_prob_thr=$(config)_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
    elseif option == "random_weight_fswp_p1_both"
        topo_params = TopoParams(
            min_dim=min_dim,
            max_dim=max_dim,
            features_list=["barcodes", "norm_barcodes", "bettis"]
        )
        # Connectivity based null model, negatve values chagned to zeros

        key_hc = "hc-both_rnd_connectivity_fswp_prob_thr$(config)"
        key_sch = "sch-both_rnd_connectivity_fswp_prob_thr$(config)"
        data_info_vector[key_hc] = DataInfo(
            key_hc,
            MatrixLoader(
                datadir("sims",
                    "COBRE_rnd_connectivity",
                ),
                "hc-both_rnd_connectivity_fswp_prob_thr=$(config)_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
        data_info_vector[key_sch] = DataInfo(
            key_sch,
            MatrixLoader(
                datadir("sims",
                    "COBRE_rnd_connectivity",
                ),
                "sch-both_rnd_connectivity_fswp_prob_thr=$(config)_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
    elseif option == "random_weight_p_all"
        topo_params = TopoParams(
            min_dim=min_dim,
            max_dim=max_dim,
            features_list=["barcodes", "norm_barcodes", "bettis"]
        )
        # Connectivity based null model, negatve values chagned to zeros
        for config in "0." .* ["05", "1", "2", "5", "9", "99"]
            key_hc = "hc_rnd_connectivity_prob_thr$(config)"
            key_sch = "sch_rnd_connectivity_prob_thr$(config)"
            data_info_vector[key_hc] = DataInfo(
                key_hc,
                MatrixLoader(
                    datadir("sims",
                        "COBRE_rnd_connectivity",
                    ),
                    "hc_rnd_connectivity_prob_thr=$(config)_",
                    "csv"),
                topo_params,
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
            )
            data_info_vector[key_sch] = DataInfo(
                key_sch,
                MatrixLoader(
                    datadir("sims",
                        "COBRE_rnd_connectivity",
                    ),
                    "sch_rnd_connectivity_prob_thr=$(config)_",
                    "csv"),
                topo_params,
                samples_limiter=samples_limiter,
                ordering_kwargs=(assign_same_values=true, ordering_start=0),
                preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
            )
        end
    elseif option == "random_weight_all_COBRE_fswp_p_x"
        topo_params = TopoParams(
            min_dim=min_dim,
            max_dim=max_dim,
            features_list=["barcodes", "norm_barcodes", "bettis"]
        )
        # Connectivity based null model, negatve values chagned to zeros


        prob_thr_str = "$(config)" |> y -> replace(y, "." => "-")
        key_pt1 = "COBRE_all_rnd_connectivity_fswp_prob_thr=$(prob_thr_str)_pt1_"
        key_pt2 = "COBRE_all_rnd_connectivity_fswp_prob_thr=$(prob_thr_str)_pt2_"
        data_info_vector[key_pt1] = DataInfo(
            key_pt1,
            MatrixLoader(
                datadir("sims",
                    "COBRE_rnd_connectivity",
                ),
                key_pt1,
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
        data_info_vector[key_pt2] = DataInfo(
            key_pt2,
            MatrixLoader(
                datadir("sims",
                    "COBRE_rnd_connectivity",
                ),
                key_pt2,
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
    elseif option == "hcp_rnd_weight_fswp_p_x"
        topo_params = TopoParams(
            min_dim=min_dim,
            max_dim=max_dim,
            features_list=["barcodes", "norm_barcodes", "bettis"]
        )
        # Connectivity based null model, negatve values chagned to zeros

        key_hcp = "hcp_rnd_connectivity_fswp_prob_thr$(config)"
        data_info_vector[key_hcp] = DataInfo(
            key_hcp,
            MatrixLoader(
                datadir("sims",
                    "HCP_rnd_connectivity",
                ),
                "HCP_1-both_rnd_connectivity_fswp_prob_thr=$(config)_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix, reverse_sign],
        )
    elseif option == "euclidean-brain_single_nx"
        noise_level = replace(data_config, "." => "-")
        topo_params = TopoParams(
            min_dim=min_dim,
            max_dim=max_dim,
            features_list=["barcodes", "norm_barcodes", "bettis"]
        )

        key_beuclidean = "brain_euclidean_rnd_n$(noise_level)"

        data_info_vector[key_beuclidean] = DataInfo(
            key_beuclidean,
            MatrixLoader(
                datadir("sims",
                    "brain-euclidean",
                    "noise=$(noise_level)"
                ),
                "brain-euclidean_rnd_connectivity_$(noise_level)_",
                "csv"),
            topo_params,
            samples_limiter=44,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix,],
        )
    elseif option == "euclidean-brain_nx"

        noise_level = replace(data_config, "." => "-")
        topo_params = TopoParams(
            min_dim=min_dim,
            max_dim=max_dim,
            features_list=["barcodes", "norm_barcodes", "bettis"]
        )


        key_beuclidean_pt1 = "brain_euclidean_pt1_rnd_n$(noise_level)"
        key_beuclidean_pt2 = "brain_euclidean_pt2_rnd_n$(noise_level)"

        data_info_vector[key_beuclidean_pt1] = DataInfo(
            key_beuclidean_pt1,
            MatrixLoader(
                datadir("sims",
                    "brain-euclidean",
                    "noise=$(noise_level)"
                ),
                "brain-euclidean_rnd_connectivity_pt1_$(noise_level)_",
                "csv"),
            topo_params,
            samples_limiter=44,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix,],
        )
        data_info_vector[key_beuclidean_pt2] = DataInfo(
            key_beuclidean_pt2,
            MatrixLoader(
                datadir("sims",
                    "brain-euclidean",
                    "noise=$(noise_level)"
                ),
                "brain-euclidean_rnd_connectivity_pt2_$(noise_level)_",
                "csv"),
            topo_params,
            samples_limiter=samples_limiter,
            ordering_kwargs=(assign_same_values=true, ordering_start=0),
            preprocessing_pipeline=Function[symmetrize_matrix],
        )
    else
        @error "Unrecognised option"
    end

    ## ===-===-
    # all_data_info = [null_model_config_, main_sch_config, main_geom_config, main_rand_config]
    # null_data_info = Dict(
    #     get_identifier(null_model_config_1) => null_model_config_1,
    #     get_identifier(null_model_config_2) => null_model_config_2,
    #     get_identifier(null_model_config_3) => null_model_config_3,
    #     get_identifier(null_model_config_4) => null_model_config_4,
    #     get_identifier(null_model_config_5) => null_model_config_5,
    #     get_identifier(null_model_config_6) => null_model_config_6,
    #     get_identifier(null_model_config_7) => null_model_config_7,
    #     get_identifier(null_model_config_8) => null_model_config_8,
    # )
    return data_info_vector
end


# end #module
