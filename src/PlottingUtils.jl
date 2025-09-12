
function yticks_formatter(values::Vector{Float64})
    max_val = max(values...)
    min_val = min(values...)
    labels = String[]
    if max_val > 5000
        labels = ["$(round(Int, value/1000))k" for value in values]
    elseif max_val > 1000
        labels = ["$(value/1000)k" for value in values]
    else
        labels = ["$(round(Int,value))" for value in values]
    end
    if min_val == 0.0
        labels[1] = "0"
    end
    return labels
end