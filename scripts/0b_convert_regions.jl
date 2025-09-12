#=
Converts default regions list into the Large scale brain networks.

=#
using DrWatson
@quickactivate "schtoppaper"

## ===-===-===-===-===-

using DataStructures: OrderedDict
using JSON
using DelimitedFiles
import Colors: RGB

## ===-===-===-===-===-
# Load 94 regions names and their coordinates
original_regions_file = "aal2_center_coords.txt"
all_regions = JSON.parsefile(datadir("exp_raw", original_regions_file)) # full path


index_to_region_name = Dict{Int,Tuple{String,Vector{Float64}}}()

all_regions_split = all_regions |> keys .|> split
for value = all_regions_split
    related_vector = all_regions[join(value, " ")]
    index_to_region_name[parse(Int, value[1])] = (value[2], related_vector)
end

total_regions = length(index_to_region_name)
## ===-===-
# Extract coordinates from regions
all_coordinates = zeros(Float64, (3, length(index_to_region_name)))

for (key, value) = index_to_region_name
    all_coordinates[:, key] .= value[2]
end

## ===-===-===-===-===-
# Load dictionary with atlases mappings
yeo7_networks_file = "AAL_Yeo_Mapping.txt"
yeo_file_content = readdlm(datadir("exp_raw", yeo7_networks_file), String)
translation_array = yeo_file_content[2:end, 1:2:3]
translation_dict = Dict()

for x in 1:94
    translation_dict[parse(Int, translation_array[x, 1])] = translation_array[x, 2]
end

for (key, value) in translation_dict
    if value == "-"
        translation_dict[key] = "Subcortical"
    end
end

## get unique large scale regions names
# This is hard coded in order to have repetability in plotting
yeo7_networks = [
    "Cont",
    "Vis",
    "Default",
    "Limbic",
    "SalVentAttn",
    "SomMot",
    "DorsAttn",
    "Subcortical"]

colours_palete =
    map(x -> RGB((x ./ 255)...),
        [
            [249, 183, 18], # Frontoparietal
            [174, 73, 177], # visual
            [233, 105, 123], # Default
            [244, 254, 195], # Limbic
            [242, 87, 255], # Ventral Attention
            [111, 155, 196],# Somatomotor
            [0, 155, 24], # Dorsal Attention
            [207, 56, 0], # Subcortical
        ]
    ) |> palette;

ls_region_to_color = OrderedDict()
for (region, color) in zip(yeo7_networks, colours_palete)
    ls_region_to_color[region] = color
end

## ===-===-===-===-===-===-===-===-
# Create a dictionary to store all information about the reigons
brain_regions_info = OrderedDict()
colors_dict = Dict()

for (key, value) in index_to_region_name
    brain_regions_info[key] = (value[1], translation_dict[key], value[2], ls_region_to_color[translation_dict[key]])
end

## ===-===-===-===-
new_regions_order = zeros(Int, total_regions)

region = yeo7_networks[1]
for region = yeo7_networks
    region_related = findall(x -> x[2] == region, brain_regions_info)
    sort!(region_related)

    first_zero_index = findfirst(x -> x == 0, new_regions_order)
    region_related_range = first_zero_index:(first_zero_index+length(region_related)-1)
    for (index, k) in enumerate(region_related_range)
        new_regions_order[k] = (region_related[index][1])
    end
end



## ===-===-===-===-===-
# Set up Yeo color scheme
regions_coloring = Array{Union{Int,Missing},2}(missing, total_regions, 1)

for (region_index, n_region) = enumerate(yeo7_networks)
    regions_related_coords = [value[3] for (key, value) in brain_regions_info if value[2] == n_region]

    large_region_relatives = [key for (key, value) in brain_regions_info if value[2] == n_region]
    regions_coloring[large_region_relatives, 1] .= region_index
end
