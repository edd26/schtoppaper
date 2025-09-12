using LsqFit # for detrending

function detrend_signal(signal::Vector, model_func::Function)
    y_vals = deepcopy(signal)
    x_vals = 1:length(y_vals)

    p0=[x_vals[1], y_vals[1]]

    fit = curve_fit(model_func, x_vals, y_vals, p0)
    y_fit = model_func(x_vals, fit.param)
    y_vals.-= y_fit

    return y_vals
end
