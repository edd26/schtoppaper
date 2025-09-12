using DrWatson
@quickactivate "SchiTopology"

## ===-===-===-
# Default configs
datainfodir(args...) = srcdir("DataInfos", args...)
"DataInfoCOBRE.jl" |> datainfodir |> include
"DataInfoHCP.jl" |> datainfodir |> include
"DataInfoGeometricModel.jl" |> datainfodir |> include
"DataInfoNull.jl" |> datainfodir |> include

# ===-===-===-
import DataStructures: OrderedDict

function configure_individuals_from_args(data_set::String, data_config::String, min_dim::Int, max_dim::Int)
    data_info = 0

    if data_set != "COBRE"
        if data_set == "geom"
            data_info = get_geometric_model_data(
                data_config;
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94
            )


        elseif data_set == "rand"
            data_info = get_geometric_model_data(
                data_config;
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94
            )
        elseif data_set == "rand"
            data_info = get_geometric_model_data(
                "1a";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94
            )
        elseif data_set == "rand_10"
            data_info = get_geometric_model_data(
                data_config;
                min_dim=0,
                max_dim=10,
                samples_limiter=44,
                matrix_size=10
            )

        elseif data_set == "Null-zero"

            data_info = get_null_model_data(
                "5";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94
            )
        elseif data_set == "Null-abs"
            data_info = get_null_model_data(
                "6";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94
            )
        elseif data_set == "Null-invgauss-zeroed"
            data_info = get_null_model_data(
                "7";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94
            )
        elseif data_set == "null_random_rewire_1"
            data_info = get_null_model_data(
                "random_rewire_1";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=100,
                matrix_size=94
            )
        elseif data_set == "null_random_rewire_2"
            data_info = get_null_model_data(
                "random_rewire_2";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=100,
                matrix_size=94
            )
        elseif data_set == "random_weight_1"
            data_info = get_null_model_data(
                "random_weight_1";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94
            )
        elseif data_set == "random_weight_p_x"
            data_info = get_null_model_data(
                "random_weight_p_x";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94,
                config=data_config
            )
        elseif data_set == "random_weight_fswp_p_x"
            data_info = get_null_model_data(
                "random_weight_fswp_p_x";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94,
                config=data_config
            )

        elseif data_set == "random_weight_fswp_p1_both"
            data_info = get_null_model_data(
                "random_weight_fswp_p1_both";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94,
                config=data_config
            )
        elseif data_set == "random_weight_all_COBRE_fswp_p_x"
            data_info = get_null_model_data(
                "random_weight_all_COBRE_fswp_p_x";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94,
                config=data_config
            )
        elseif data_set == "COBRE_vs_random_weight_all_COBRE_fswp_p_x"
            data_info_rand_all = get_null_model_data(
                "random_weight_all_COBRE_fswp_p_x";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94,
                config=data_config
            )
            rand_keys = data_info_rand_all |> keys |> collect

            data_info = OrderedDict(
                rand_keys[1] => data_info_rand_all[rand_keys[1]],
                rand_keys[2] => data_info_rand_all[rand_keys[2]],
                "hc" => get_COBRE_data_info("hc"),
                "sch" => get_COBRE_data_info("sch"),
            )

        elseif data_set == "random_weight_p_all"
            data_info = get_null_model_data(
                "random_weight_p_all";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94,
                config=data_config
            )


        elseif data_set == "HCP_hc"
            data_info = OrderedDict(
                "HCP_test_orig" => get_HCP_config("HCP_test_orig"),
                "hc" => get_COBRE_data_info("hc"),
            )
        elseif data_set == "HCP1_hc"
            data_info = OrderedDict(
                "HCP_1" => get_HCP_config("HCP_1", samples_limiter=44),
                "hc" => get_COBRE_data_info("hc"),
            )
        elseif data_set == "HCP1_sch"
            data_info = OrderedDict(
                "HCP_1" => get_HCP_config("HCP_1", samples_limiter=44),
                "sch" => get_COBRE_data_info("sch"),
            )
        elseif data_set == "COBRE-not-ordered"
            data_info = OrderedDict(
                "hc_not_orded" => get_COBRE_data_info("hc_not_orded"),
                "sch_not_orded" => get_COBRE_data_info("sch_not_orded"),
            )
        elseif data_set == "HCP_sch"
            data_info = OrderedDict(
                "HCP_test_orig" => get_HCP_config("HCP_test_orig"),
                "sch" => get_COBRE_data_info("sch"),
            )
        elseif data_set == "HCP_1"
            data_info = OrderedDict(
                "HCP_1" => get_HCP_config("HCP_1"),
            )
        elseif data_set == "HCP_1_rev"
            data_info = OrderedDict(
                "HCP_1_rev" => get_HCP_config("HCP_1_rev"; max_dim=max_dim, min_dim=min_dim),
            )
        elseif data_set == "HCP_1_modified"
            data_info = OrderedDict(
                "hc" => get_COBRE_data_info("hc"),
                "sch" => get_COBRE_data_info("sch"),
                "HCP_1_modified" => get_HCP_config("HCP_1_modified"),
            )
        elseif data_set == "COBRE_HCP"
            data_info = OrderedDict(
                "hc" => get_COBRE_data_info("hc"),
                "sch" => get_COBRE_data_info("sch"),
                "HCP_1" => get_HCP_config("HCP_1"),
            )

        elseif data_set == "hcp_rnd_weight_fswp_p_x"
            data_info = get_null_model_data(
                "hcp_rnd_weight_fswp_p_x";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=99,
                matrix_size=94,
                config=data_config
            )
        elseif data_set == "euclidean-brain_single_nx"
            data_info = get_null_model_data(
                "euclidean-brain_single_nx";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94,
                config=data_config
            )
        elseif data_set == "euclidean-brain_nx"
            data_info = get_null_model_data(
                "euclidean-brain_nx";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94,
                config=data_config
            )

        elseif data_set == "COBRE_random_weight_fswp_p_x"
            data_info_cobre = [
                "hc" => get_COBRE_data_info("hc"),
                "sch" => get_COBRE_data_info("sch"),
            ] |> OrderedDict

            data_info_fswp = get_null_model_data(
                "random_weight_fswp_p_x";
                min_dim=min_dim,
                max_dim=max_dim,
                samples_limiter=44,
                matrix_size=94,
                config=data_config
            )
            data_info = vcat([
                [k => v for (k, v) in data_info_fswp]...,
                [k => v for (k, v) in data_info_cobre]...,
            ]) |> OrderedDict
        elseif data_set == "COBRE_HCP_random_weight_fswp_p_all"
            data_info_cobre = [
                "hc" => get_COBRE_data_info("hc"),
                "sch" => get_COBRE_data_info("sch"),
            ] |> OrderedDict

            all_configs = ["$(c)" for c in [0.05, 0.1, 0.2, 0.5, 0.9, 0.99]]
            all_data_info_fswp_pt1 = [
                get_null_model_data(
                    "random_weight_fswp_p_x";
                    min_dim=min_dim,
                    max_dim=max_dim,
                    samples_limiter=44,
                    matrix_size=94,
                    config=c
                ) for c in all_configs]
            all_data_info_fswp = vcat(
                [vcat([v for v in s]...) for s in all_data_info_fswp_pt1]...
            ) |> OrderedDict

            data_info = vcat([
                [k => v for (k, v) in data_info_cobre]...,
                [k => v for (k, v) in all_data_info_fswp]...,
            ]) |> OrderedDict
            data_info["HCP_1"] = get_HCP_config("HCP_1")

        elseif occursin("HCP", data_set) && occursin("1", data_set)
            data_info = OrderedDict(
                data_set => get_HCP_config(data_set),
            )
        elseif occursin("HCP", data_set)
            data_info = OrderedDict(
                data_set => get_HCP_config(data_set),
            )

        elseif data_set == "COBRE_rev"
            data_info = OrderedDict(
                "rev_hc" => get_COBRE_data_info("rev_hc", min_dim=min_dim, max_dim=max_dim),
                "rev_sch" => get_COBRE_data_info("rev_sch", min_dim=min_dim, max_dim=max_dim),
            )
        elseif data_set == "COBRE_patch_weak"
            data_info = OrderedDict(
                "patched_weak_hc" => get_COBRE_data_info("patched_weak_hc", min_dim=min_dim, max_dim=max_dim),
                "patched_weak_sch" => get_COBRE_data_info("patched_weak_sch", min_dim=min_dim, max_dim=max_dim),
            )
        else
            "Did not find any dataset that matched arg param" |> ErrorException |> throw
        end

    else
        data_info = OrderedDict(
            "hc" => get_COBRE_data_info("hc", min_dim=min_dim, max_dim=max_dim),
            "sch" => get_COBRE_data_info("sch", min_dim=min_dim, max_dim=max_dim),
        )

    end # if

    return data_info
end
