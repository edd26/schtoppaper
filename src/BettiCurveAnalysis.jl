
function interpolate_and_average(vectors)
    # Find the maximum length among the vectors
    max_length = maximum(length.(vectors))

    # Create a new range from 0 to 1 with max_length points
    x_new = range(0, stop=1, length=max(2, 2 * max_length))

    # Initialize arrays to store interpolated vectors
    interpolated_vectors = [zeros(max_length) for _ in 1:length(vectors)]

    # Iterate through each input vector
    for (i, v) in enumerate(vectors)
        # Create interpolation object
        itp = LinearInterpolation(range(0, stop=1, length=max(2, length(v))), v)

        # Evaluate the interpolation at the new x-values
        interpolated_vectors[i] = itp(x_new)
    end

    # Calculate the average
    avg_vector = sum(interpolated_vectors) ./ length(vectors)

    return interpolated_vectors, avg_vector
end


@enum RegionSelection begin
    Left
    Right
    Both
end


function get_all_required_data_structures(rand_matrices_info, rand_symmetric_matrices; min_dim::Int=0, max_dim::Int=3)
    ordering_kwargs_rand = rand_matrices_info |> get_ordering_kwargs
    total_matrices = size(rand_symmetric_matrices, 1)

    ordered_rand_matrices = [TopologyPreprocessing.get_ordered_matrix(
        rand_symmetric_matrices[k, :, :];
        ordering_kwargs_rand...)
                             for k in 1:total_matrices]
    rand_matrices_C =
        [C = eirene(mat, maxdim=max_dim, model="vr")
         for mat in ordered_rand_matrices]

    rand_params = [
        get_topo_features(C, topo_params |> get_features, min_dim=min_dim, max_dim=max_dim)
        for C in rand_matrices_C]

    return ordered_rand_matrices, rand_matrices_C, rand_params
end
