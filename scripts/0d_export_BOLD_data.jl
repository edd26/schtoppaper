#=
Load original numpy files with BOLD signal and connectivity matrices and exports
all matrices as separate CSV files.

=#

module BOLDdata


export FC_dir, T_dir
## ===-===-
import DrWatson: @quickactivate, datadir, srcdir
@quickactivate "schtoppaper"

## ===-===-
using NPZ # for npy files
using Pipe
using LsqFit # for detrending

using DelimitedFiles
using Statistics

"preprocessing_utils.jl" |> srcdir |> include
## ===-===-
data_keys = ["hc", "sch"]
bold_keys = data_keys .* "_RSC"

FC_dir(args...) = datadir("exp_pro", "AvgFmatrix", args...)
T_dir(args...) = datadir("exp_pro", "AvgTmatrix", args...)

functional_files_location = @pipe ["AvgFmatrixHC.npy"; "AvgFmatrixSCZ.npy"] .|>
                                  datadir("exp_raw", _)

time_series_files_location = @pipe ["AvgTmatrixHC.npy"; "AvgTmatrixSCZ.npy"] .|>
                                   datadir("exp_raw", _)

## ===-===-
average_FC_matrix = Dict()
for (path, key) in zip(functional_files_location, bold_keys)
    average_FC_matrix[key] = npzread(path)
end

average_T_matrix = Dict()
for (path, key) in zip(time_series_files_location, bold_keys)
    average_T_matrix[key] = npzread(path)
end


## ===-
# Do the detrending for time series
detrended_signals = deepcopy(average_T_matrix)
lin_model(x, a) = a[1] .* x .+ a[2]

for key in bold_keys
    total_subjects, total_regions, total_samples = size(average_T_matrix[key])

    for subject = 1:total_subjects, region = 1:total_regions
        detrended_signals[key][subject, region, :] = detrend_signal(
            average_T_matrix[key][subject, region, :],
            lin_model
        )
    end
end

## ===-
# Compute pearson correlation between regions within a subject
detrended_pcorr = deepcopy(average_FC_matrix)

for key in bold_keys
    total_subjects, total_regions, total_samples = size(average_T_matrix[key])

    for subject = 1:total_subjects
        y1_vals = detrended_signals[key][subject, :, :]
        y2_vals = detrended_signals[key][subject, :, :]

        detrended_pcorr[key][subject, :, :] = cor(y1_vals, y2_vals, dims=2)
    end
end

end # module
