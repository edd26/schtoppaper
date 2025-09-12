srcdir("DataSource", "SignalSource.jl") |> include

## ===-===-
"""
    SignalGenerator <: SignalSource

Information about how to generate matrices used as a sequence

"""
struct SignalGenerator <: SignalSource

    samples_count::Int
    generation_function::Function
    generation_arguments::Tuple

    matrix_size::Int
    total_matrices::Int

    function SignalGenerator(samples_count, generation_function, generations_arguments; matrix_size=94, total_matrices=104)
        new(samples_count, generation_function, generations_arguments, matrix_size, total_matrices)
    end

end

get_samples_count(signal_matrix::SignalGenerator) = signal_matrix.samples_count
get_generator_function(signal_matrix::SignalGenerator) = signal_matrix.generation_functio
get_generator_arguments(signal_matrix::SignalGenerator) = signal_matrix.generation_arguments

"""

This function return a structure with only shape of the desired data; the data itself
    is generated when fragment correlation should be computed.

# TODO this should be refactored so that this information is taken from the structure
    and not empty data
"""
function get_matrices(sig_mat_info::SignalGenerator;)
    func = sig_mat_info.generation_function
    args = sig_mat_info.generation_arguments

    total_sessions = sig_mat_info.samples_count
    matrix_size = sig_mat_info.matrix_size
    total_matrices = sig_mat_info.total_matrices

    final_matrix = zeros(total_sessions, matrix_size, total_matrices)

    return final_matrix
end


function get_matrix_identifier(sig_mat_info::SignalGenerator)
    id = "samples=$(mat_gen.samples_count)_matrix_size=$(mat_gen.matrix_size)"
end
