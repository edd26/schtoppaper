module SchiArgPar

using ArgParse

function parse_clustering_commandline()
    s = ArgParseSettings()

    # ===
    @add_arg_table! s begin

        "--data_set", "-d"
        help = "Specifies data set to use"
        arg_type = String
        default = "COBRE" # COBRE
        # default = "COBRE_rev" # COBRE
        # default = "HCP_1"
        # default = "HCP_1_modified"
        # default = "COBRE_HCP"
        # default = "HCP_1_rev"
        # default = "geom"
        # default = "rand"
        # default = "rand_10"
        # default = "random_weight_1"
        # default = "random_weight_p_x"
        # default = "random_weight_p_all"
        # default = "random_weight_fswp_p_x"
        # default = "random_weight_all_COBRE_fswp_p_x"
        # default = "COBRE_vs_random_weight_all_COBRE_fswp_p_x"
        # default = "hcp_rnd_weight_fswp_p_x"
        # default = "random_weight_fswp_p1_both"
        # default = "euclidean-brain_nx"
        # default = "euclidean-brain_single_nx"
        # default = "COBRE_random_weight_fswp_p_x"
        # default = "COBRE_HCP_random_weight_fswp_p_all"


        "--data_config", "-c"
        help = "Specifies data config to use (see src/DataInfo* files for possible options)"
        arg_type = String
        default = "1e"

        "--min_dim"
        help = "Minimal dimension for topological computations"
        arg_type = Int
        default = 0

        "--max_dim"
        help = "Maximal dimension for topological computations"
        arg_type = Int
        default = 3

        "--selected_dim", "-s"
        help = "Dimension selected for clustering analysis "
        arg_type = Int
        default = 1

        "--limit_popularity"
        help = "Sets lower limit of the cycle popularity "
        arg_type = Int
        default = 0

        "--total_shuffles"
        help = "Total number of shuffles for distribution generation"
        arg_type = Int
        default = 10000


        "--force_symmat_loading"
        help = "Enables overwriting of files with symmetrised matrices."
        arg_type = Bool
        default = false

        "--force_ordering"
        help = "Enables overwriting of files with ordered matrices."
        arg_type = Bool
        default = false

        "--force_topology"
        help = "Enables overwriting of files with topological data."
        arg_type = Bool
        default = false

        "--force_topo_features"
        help = "Enables overwriting of files with topological data."
        arg_type = Bool
        default = false

        "--force_computations"
        help = "Enables overwriting of files in produce_or_load."
        arg_type = Bool
        default = false

        "--force_landscpes_permutation"
        help = "Enables overwriting of landscapes permutation in 8e script"
        arg_type = Bool
        default = false

        "--force_clusters_dist_permutation"
        help = "Enables overwriting of distance distribution computations in 6cblcb"
        arg_type = Bool
        default = false


        "--normalisation_type"
        help = "Decides on the type of normalisation, also affecting orderig type. 3 options are possible: 0 (default) is for normalising to number of unique values (ordering is done without step); 1 is sports normalising, ordering is done with step and normalisation is done with last placement; 2 is max normalisation with total number of entries in upper diagonal of matrix (ordering is done with a step)"
        arg_type = Int
        default = 0

        "--use_outer_layer_only"
        help = "For landscapes operations, decides if use outer layer only or all"
        arg_type = Bool
        default = false

        "--linkage", "-l"
        help = "Specifies linkage to use"
        arg_type = String
        default = "ward"
        # required = true

        "--min_clusters", "-t"
        help = "Specifies the maximal height for the clusters"
        arg_type = Int
        default = 3

        "--cluster_height", "-w"
        help = "Specifies the maximal height for the clusters"
        arg_type = Int
        default = 0

        "--p_value", "-p"
        help = "p value for p-norm computations"
        arg_type = Int
        default = 1

        "--zero_landscapes"
        help = "Use landscapes that are only zero for distribution with less than 44 subjects per group"
        arg_type = Bool
        default = true

        "--wasserstein_dsitance"
        help = "Bool condition for using the Wasserstein metri, default false"
        arg_type = Bool
        default = false

        "--manual_clustering"
        help = "If different from 0, then applies n-split of the cluster, starting from the top"
        arg_type = Int
        default = 0

        "--minimal_height"
        help = "If manual clustering is not 0, then this specifies the minimal height for the clusters"
        arg_type = Int
        default = 0

        "--minimal_popularity"
        help = "If manual clustering is not 0, then this specifies the minimal number of members in the cluster"
        arg_type = Int
        default = 0

        "--allow_final_splits"
        help = "If manual clustering is not 0, then this specifies the minimal number of members in the cluster"
        arg_type = Bool
        default = true

        "--max_width"
        help = "If manual clustering is not 0, then this specifies the maximal number of cycles in a cluster"
        arg_type = Int
        default = 300

        "--min_width"
        help = "If manual clustering is not 0, then this specifies the minimal number of cycles in a cluster"
        arg_type = Int
        default = 20

        "--do_vector_extension"
        help = "If set true, extends subjects vectors with information about which Yeo network is present in the vector."
        arg_type = Bool
        default = false

        "--repeat_cycles"
        help = "If set true, repeats the cycle as many times as it is popular (was default behaviour before introduced)."
        arg_type = Bool
        default = false

        "--show_signifficant"
        help = "If set true, produces results for the cycles that are found significantly different between the groups. Used in 6cblcbbgh tex-export script."
        arg_type = Bool
        default = true

    end

    return parse_args(s)
end
end # module

#
# Data set export:

function get_data_info_from_input(data_set, data_config, samples_limiter)
    all_sets = split(data_set, ",") .|> String
    all_configs = split(data_config, ",") .|> String
    all_limiters = split(samples_limiter, ",") .|> String
    return ["$(k)-$(v)" => Dict("set" => k, "conf" => v, "limiter" => l) for (k, v, l) in zip(all_sets, all_configs, all_limiters)] |> OrderedDict
end


function get_data_info(data_set, data_config, min_dim, max_dim; samp_lim="1", shuffling::Symbol=:full)
    samples_limiter = parse(Int, samp_lim)
    if data_set == "BOLD"
        data_info = get_BOLD_data_set(data_config;
            min_dim=min_dim,
            max_dim=max_dim,
            samples_limiter=samples_limiter
        )
    elseif data_set == "CRT0"
        data_info = get_CRT0_data_set(data_config;
            min_dim=min_dim,
            max_dim=max_dim
        )
    elseif data_set == "COBRE"
        data_info = get_COBRE_data_set(data_config;
            min_dim=min_dim,
            max_dim=max_dim
        )
    elseif data_set == "COBRE_shuff"
        data_info = get_COBRE_shuffled_data_set(data_config;
            min_dim=min_dim,
            max_dim=max_dim,
            samples_limiter=samples_limiter,
            shuffling=shuffling
        )
    elseif data_set == "COBRE_tshuff"
        data_info = get_COBRE_time_shuffled_data_set(data_config;
            min_dim=min_dim,
            max_dim=max_dim,
            samples_limiter=samples_limiter
        )
    elseif data_set == "Geom"
        data_info = get_geometric_data_set(data_config;
            min_dim=min_dim,
            max_dim=max_dim,
            samples_limiter=samples_limiter,
            matrix_size=total_nodes,
            simulation_len=simulation_len
        )
    elseif data_set == "HRF"
        data_info = get_hrf_data_set(data_config;
            min_dim=min_dim,
            max_dim=max_dim,
            samples_limiter=samples_limiter
            # matrix_size=total_nodes
            # simulation_len=simulation_len
        )
    else
        @warn "Unrecognised data set"
        data_info = Dict()
    end
    return data_info
end
