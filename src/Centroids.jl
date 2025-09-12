
function get_landscape_centroid(pland)
    # Corner case
    if length(pland.land[1]) == 1
        x_vals, y_vals, x_hat, y_hat = 0, 0, 0, 0
        return x_vals, y_vals, x_hat, y_hat
    end
    x_vals = [l.first for l in pland.land[1]]
    y_vals = [l.second for l in pland.land[1]]
    x_hat = sum(x_vals .* y_vals) / sum(y_vals)
    y_hat = sum(y_vals .^ 2) / sum(2 * y_vals)
    return x_vals, y_vals, x_hat, y_hat
end

function cart_to_polar(x::Float64, y::Float64)
    r = sqrt(x^2 + y^2)
    θ = atan(y, x)

    return r, θ
end

function get_arrow_colour(arrow_head_x, arrow_head_y)
    if arrow_head_x[1] <= 0 && arrow_head_y[1] <= 0
        return :red
    elseif arrow_head_x[1] >= 0 && arrow_head_y[1] >= 0
        return :green
    elseif arrow_head_x[1] > 0 && arrow_head_y[1] < 0
        return :black
    elseif arrow_head_x[1] < 0 && arrow_head_y[1] > 0
        return :gray
    else
        ErrorException("Out of scope for arrow colour assigment") |> throw
    end
end