import DrWatson: srcdir
@quickactivate "SchiTopology"
srcdir("DataSource", "SymmetricMatrixSource.jl") |> include

"""
    MatrixLoader <: SymmetricMatrixSource

Information about how to load symmetric matrices.

Upon creation of this structure:
    - 'data_path' is checked for existence;
    - 'data_path'/'name_signature'k.'file_extension' is checked for existenece, where k is
    file index and k=1 for this check
    - 'samples_count' is set to max(k) of 'data_path'/'name_signature'k.'file_extension'
    - 'name_signature' must not end with a number, because this will fail the assertion

It is assumed that files indices are located at the end of the 'name_signature'.

"""
struct MatrixLoader <: SymmetricMatrixSource
    data_path::String
    name_signature::String
    file_extension::String

    samples_count::Int
    samples_iterator::Vector{String}

    function MatrixLoader(data_path::String, name_signature::String, file_extension::String; samples_count::Int=0)
        # if there is a number in signiture, then loading will fail
        check_number_ending(name_signature)

        fname = .*(name_signature, "1.", file_extension)
        first_file = joinpath(data_path, fname)

        # To ensure correct loading:
        check_if_exsits(data_path, ispath)
        check_if_exsits(first_file, isfile)

        # To have correct samples_count
        samples_iterator = get_sample_iterator(data_path, name_signature, file_extension)
        if samples_count == 0
            samples_count = length(samples_iterator)
        end

        new(data_path, name_signature, file_extension, samples_count, samples_iterator)
    end
end

get_path(mat_gen::MatrixLoader) = mat_gen.data_path
get_signature(mat_gen::MatrixLoader) = mat_gen.name_signature
get_extension(mat_gen::MatrixLoader) = mat_gen.file_extension
get_samples_count(mat_gen::MatrixLoader) = mat_gen.samples_count
get_iterator(mat_gen::MatrixLoader) = mat_gen.samples_iterator

"""
Creates a dictionary with arguments used for file generaion.
"""
function get_matrix_identifier(mat_load::MatrixLoader)
    return mat_load.name_signature
end

function get_matrices(matrix_info::MatrixLoader; default_read_args=(',', Float64, '\n'), counter_start::Int=1)
    fpath = get_path(matrix_info)
    fprefix = get_signature(matrix_info)
    fextension = get_extension(matrix_info)
    total_samples = get_samples_count(matrix_info)

    get_full_name(file_name) = joinpath(fpath, file_name)
    read_data(fname) = readdlm(fname, default_read_args...)

    all_matrices = Any[]
    if counter_start == 0
        index_range = counter_start:(total_samples-1)
    else
        index_range = counter_start:total_samples
    end


    for index = index_range
        push!(all_matrices, [fprefix, "$(index).", fextension] |> join |> get_full_name |> read_data)
    end
    total_regions, total_points = size(all_matrices[1])
    final_matrix = zeros(total_samples, total_regions, total_points)
    for (position, index) = index_range |> enumerate
        final_matrix[position, :, :] = all_matrices[position]
    end
    return final_matrix
end



function get_matrix(sym_mat_info::MatrixLoader, index::Int; default_read_args=(',', Float64, '\n'))
    fpath = get_path(sym_mat_info)
    fprefix = get_signature(sym_mat_info)
    file_extension = get_extension(sym_mat_info)

    get_full_name(file_name) = joinpath(fpath, file_name)
    read_data(fname) = readdlm(fname, default_read_args...)

    return [fprefix, "$(index).", file_extension] |> join |> get_full_name |> read_data
end
