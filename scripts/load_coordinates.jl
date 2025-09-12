"helper_functions.jl" |> srcdir |> include

## =====-=====-=====-=====-=====-=====-=====-=====-=====-=====-=====-=====-=====-
# Get number-region dictionary
if !@isdefined file_name
    file_name = "aal2_center_coords.txt"
    number_to_region = get_regions_from_file(file_name)

    region_coordinates = get_region_coords(file_name)
else
    @info "Using file name form workspace"
end

# Plot brain with regions
## ===-===-
# Create matrix with all coordinates
if !@isdefined all_coordinates
    all_coordinates = zeros(Float64, (3, length(region_coordinates)))

    for (index, (region, coordinates)) in enumerate(region_coordinates)
        all_coordinates[:, index] .= coordinates
    end
else
    @info "Using all_coordinates form workspace"
end

if !@isdefined region_label
    region_label = String[]

    for (index, (region, coordinates)) in enumerate(region_coordinates)
        push!(region_label, region)
    end
else
    @info "Using region_label form workspace"
end
