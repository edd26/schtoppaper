import DrWatson: srcdir
@quickactivate "SchiTopology"
srcdir("DataSource", "SymmetricMatrixSource.jl") |> include


"""
    MatrixGenerator <: SymmetricMatrixSource(...

Stores information about how to run 'generator_function' 'samples_count' many
times with 'generator_args'.

'samples_count' must be a positive integer
"""
##
struct MatrixGenerator <: SymmetricMatrixSource
    samples_count::Int
    matrix_size::Int
    generator_function::Function

    generator_args::Union{Dict,Tuple}
    generator_kwargs::Union{Dict,NamedTuple}

    function MatrixGenerator(
        samples_count::Int,
        matrix_size::Int,
        generator_function::Function;
        generator_args::Tuple=(NaN,),
        generator_kwargs::NamedTuple=NamedTuple()
    )

        new(
            samples_count,
            matrix_size,
            generator_function,
            generator_args,
            generator_kwargs,
        )
    end

    function MatrixGenerator(
        samples_count::Int,
        matrix_size::Int,
        generator_function::Function,
        generator_args::Union{Dict,Tuple},
        generator_kwargs::Union{Dict,NamedTuple},
    )
        samples_count == 0 && DomainError(samples_count, "Can not run generator 0 times.") |> throw
        samples_count <= 0 && DomainError(samples_count, "Can not run generator for negative number.") |> throw

        matrix_size == 0 && DomainError(matrix_size, "Can not run generator for matrix of size 0.") |> throw
        matrix_size <= 0 && DomainError(matrix_size, "Can not run generator for matrix size being negative number.") |> throw

        new(
            samples_count,
            matrix_size,
            generator_function,
            generator_args,
            generator_kwargs
        )
    end
end
##
get_function(mat_gen::MatrixGenerator) = mat_gen.generator_function
get_args(mat_gen::MatrixGenerator) = mat_gen.generator_args
get_kwargs(mat_gen::MatrixGenerator) = mat_gen.generator_kwargs
get_samples_count(mat_gen::MatrixGenerator) = mat_gen.samples_count

"""
Creates a dictionary with arguments used for file generaion.
"""
function get_matrix_identifier(mat_gen::MatrixGenerator)
    id = "samples=$(mat_gen.samples_count)_matrix_size=$(mat_gen.matrix_size)"

    return @strdict id
end


## ===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-
function get_symmetric_matrix(container::Dict)
    @unpack data_info, index = container

    loaded_matrix = @pipe data_info |>
                          get_matrix_source |>
                          get_matrix(_, index)

    # check symmetry
    if !issymmetric(loaded_matrix)
        @info "Loaded matrix is not symmetric"
    end

    # preprocess
    preprocessing_functions = get_preprocessing(data_info)
    if !isempty(preprocessing_functions)
        @info "Preprocessing..."
        processed_matrix = loaded_matrix
        for func in preprocessing_functions
            processed_matrix = processed_matrix |> func
        end
        final_matrix = processed_matrix
    else
        final_matrix = loaded_matrix
    end

    symmetry_check = final_matrix |> !issymmetric
    symmetry_check && ("Matrix is not symmetric after processin! Check the pipeline. Aborting." |> AssertionError |> throw)

    return Dict("symmetric_matrix" => final_matrix)
end

function get_matrix(sym_mat_info::MatrixGenerator, index::Int)
    func = get_function(sym_mat_info)
    func_args = get_args(sym_mat_info)

    return func(func_args...)
end

function get_matrix(sym_mat_info::MatrixGenerator)
    func = get_function(sym_mat_info)
    func_args = get_args(sym_mat_info)

    return func(func_args...)
end

function get_matrices(sig_mat_info::MatrixGenerator;)
    total_matrices = sig_mat_info.samples_count
    matrix_size = sig_mat_info.matrix_size
    func = sig_mat_info.generator_function
    if isnan(sig_mat_info.generator_args[1])
        args = (matrix_size,)
    else
        args = sig_mat_info.generator_args
        args = (matrix_size, args...,)
    end


    matrices_vec = [func(args...) for k in 1:total_matrices]
    final_matrix = zeros(total_matrices, matrix_size, matrix_size)
    for k in 1:total_matrices
        final_matrix[k, :, :] = matrices_vec[k]
    end

    return final_matrix
end

