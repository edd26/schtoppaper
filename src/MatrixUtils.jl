

function inverse_elements(matrix_set)
    inversed_matrix = deepcopy(matrix_set)
    for (index, a) in enumerate(matrix_set)
        inversed_a = 1 ./ a
        all_infs = findall(x -> x == Inf || x == -Inf, inversed_a)
        if !isempty(all_infs)
            @debug "Replaced Infinity elements with zeros"
            inversed_a[all_infs] .= 0
            inversed_matrix[index] = inversed_a
        end
    end
    return inversed_matrix
end

# ===-===-===-===-===-
function normalise_matrix(matrix_set)
    normed_matrix = deepcopy(matrix_set)
    for (index, mat) in enumerate(matrix_set)
        min_val = findmin(mat)[1]
        max_val = findmax(mat)[1]

        normed_matrix[index] = (mat .- min_val) ./ (max_val .- min_val)
    end
    return normed_matrix
end

function standarize_matrix(mat)
    mean_val = mean(mat)[1]
    std_val = std(mat)

    standarized_matrix = (mat .- mean_val) ./ (std_val)

    return standarized_matrix
end
