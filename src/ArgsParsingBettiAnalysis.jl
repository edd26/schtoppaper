module ArgParseBettis

using ArgParse

function parse_betti_commandline()
    s = ArgParseSettings()

    # ===
    @add_arg_table! s begin

        "--min_dim"
        help = "Minimal dimension for topological computations"
        arg_type = Int
        default = 0

        "--max_dim"
        help = "Maximal dimension for topological computations"
        arg_type = Int
        default = 2

        "--samples_limiter"
        help = "Total matrices per dataset"
        arg_type = Int
        default = 44

        "--matrix_size"
        help = "Size of a single matrix"
        arg_type = Int
        default = 24

        "--regions_side"
        help = "Selects which indices from COBRE and HCP are used: 0- left regions, 1- right regions, 2- both sides"
        arg_type = Int
        default = 2

        "--highest_shuffling"
        help = "Set the uppee limit for shuffling of the matrix, expressed in a fraction from 0 to 1"
        arg_type = Float64
        default = 0.4

        "--work_on_reverse"
        help = "Controls if the computations should be done for forward or reverse filtration. Warning: reverse filtration is more computationaly expensive."
        arg_type = Bool
        default = true
    end

    return parse_args(s)
end
end # module
