srcdir("DataSource", "SignalSource.jl") |> include

using Pipe

## ===-===-
"""
    SignalMatrixLoader <: SignalSource

Information about how to load matrix with signals.

Upon creation of this structure:
    - 'data_path' is checked for existence;
    - 'data_path'/'name_signature'k.'file_extension' is checked for existenece, where k is
    file index and k=1 for this check
    - 'samples_count' is set to max(k) of 'data_path'/'name_signature'k.'file_extension'
    - 'name_signature' must not end with a number, because this will fail the assertion

It is assumed that files indices are located at the end of the 'name_signature'.

"""
struct SignalMatrixLoader <: SignalSource
    data_path::String
    name_signature::String
    file_extension::String

    samples_count::Int
    samples_iterator::Vector{String}

    function SignalMatrixLoader(data_path::String, name_signature::String, file_extension::String)
        # if there is a number in signiture, then loading will fail
        check_number_ending(name_signature)

        fname = .*(name_signature, "1.", file_extension)
        first_file = joinpath(data_path, fname)

        # To ensure correct loading:
        check_if_exsits(data_path, ispath)
        check_if_exsits(first_file, isfile)

        samples_iterator = get_sample_iterator(data_path, name_signature, file_extension)
        samples_count = length(samples_iterator)

        new(data_path, name_signature, file_extension, samples_count, samples_iterator)
    end
end

get_path(signal_matrix::SignalMatrixLoader) = signal_matrix.data_path
get_signature(signal_matrix::SignalMatrixLoader) = signal_matrix.name_signature
get_extension(signal_matrix::SignalMatrixLoader) = signal_matrix.file_extension
get_samples_count(signal_matrix::SignalMatrixLoader) = signal_matrix.samples_count
get_iterator(signal_matrix::SignalMatrixLoader) = signal_matrix.samples_iterator

function get_matrices(sig_mat_info::SignalMatrixLoader; default_read_args=(',', Float64, '\n'))
    fpath = get_path(sig_mat_info)
    fprefix = get_signature(sig_mat_info)
    fextension = get_extension(sig_mat_info)
    total_samples = get_samples_count(sig_mat_info)

    get_full_name(file_name) = joinpath(fpath, file_name)
    read_data(fname) = readdlm(fname, default_read_args...)

    all_matrices = Any[]
    for index = 1:total_samples
        push!(all_matrices, [fprefix, "$(index).", fextension] |> join |> get_full_name |> read_data)
    end
    total_regions, total_points = size(all_matrices[1])
    final_matrix = zeros(total_samples, total_regions, total_points)
    for index = 1:total_samples
        final_matrix[index, :, :] = all_matrices[index]
    end
    return final_matrix
end
