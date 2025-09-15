# Mesoscale differences in brain organization


The code provided here reproduces the figures from the paper:

* E. Dmitruk, Ch. Metzner, V. Steuber, S. N. Kadir: [Mesoscale differences in brain organization in schizophrenia revealed by topological data analysis](https://doi.org/10.1101/2025.06.19.660631) (p. 2025.06.19.660631v1). bioRxiv.

## Introduction
The code was written by Emil Dmitruk*, and the underlying ideas are the result of joint work with [Shabnam Kadir(*)](https://github.com/shabnamkadir),  Christoph Metzner and Volker Steuber(*).

(*)[UHBiocomputation Group](http://biocomputation.herts.ac.uk/)

## Code
This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> schtoppaper

It is authored by Emil Dmitruk.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "schtoppaper"
```
which auto-activate the project and enable local path handling from DrWatson.

## Data

Processed data can be obtained from our [GIN repository](https://gin.g-node.org/dreamy1494/schtoppaper-data).

# Reproducing results

*Results were reproduced on MacOS (Intel, ARM), Linux. Windows system was not tested*

## Running scripts

The scripts were designed to be run from the command line. Below we present how they can be used.
All arguments, with short descriptions are outlined in `src/ArgsParsing*.jl` files.

### Reproducing figures

In order to replicate figures, the following commands have to be run:

1. Figure 1:
```julia
julia scripts/1ga2_plot_landcscapes_samples_paper_resized.jl
```
2. Figure 2:
```julia
julia scripts/6cblcba7_Makie_plot_cluster_structure_individual_ppl_popularity_vbarcodes.jl -d COBRE --selected_dim 0 --popularity_limit 0
julia scripts/6cblcba7_Makie_plot_cluster_structure_individual_ppl_popularity_vbarcodes.jl -d COBRE --selected_dim 1 --popularity_limit 0
```
3. Figure 3:
```julia
julia scripts/8b3_makie_plot_average_persistence_landscapes_resized.jl -d COBRE --selected_dim 0 --popularity_limit 0
```
4. Figure 4:
```julia
julia scripts/6cblcbc4a_Makie_plot_landscapes_centroids_same_canvas_with_clustering.jl -d COBRE --selected_dim 1 --popularity_limit 22
```
5. Figure 5:
```julia
julia scripts/6cblcba7_Makie_plot_cluster_structure_individual_ppl_popularity_vbarcodes.jl -d COBRE_HCP --selected_dim 0 --popularity_limit 0
julia scripts/6cblcba7_Makie_plot_cluster_structure_individual_ppl_popularity_vbarcodes.jl -d COBRE_HCP --selected_dim 1 --popularity_limit 0
```
6. Figure 6:
```julia
julia scripts/6cb_get_cycles_posets.jl -d COBRE --selected_dim 0 --popularity_limit 0 --max_dim 2
julia scripts/6cb_get_cycles_posets.jl -d COBRE_rev --selected_dim 0 --popularity_limit 0 --max_dim 2
julia scripts/6cb_get_cycles_posets.jl -d HCP --selected_dim 1 --popularity_limit 0 --max_dim 2
julia scripts/6cb_get_cycles_posets.jl -d HCP_rev  --selected_dim 0 --popularity_limit 0 --max_dim 2
julia scripts/2j2_small_matrix_Betti_curves_investigation.jl --max_dim 2
```
- please note that maximal dimension is limited to dimension 2; otherwise, significant resources (RAM memory, disk space) are required to compute results
presented in the paper (up to dimension 4)
7. Figure 7:
```julia
julia scripts/1k2_matrix_comparison_plts_with_symmetric.jl
```
8. Figure 8:
```julia
julia scripts/1gb3-2_plot_landcscape_evolution_resized_paper_submission.jl
```

## Containers

We provide a `Dockerfile.julia1-8` that ensures the results can be reproduced. Below 
are the instructions how to use 2 popular contarization tools.

### Podman

1. Build:
```Podman
  podman build -t schtoppaper:latest -f Dockerfile.julia1-8
```
- `-t` is used to name the container
- `f` indicates the Docker file to use in image creation

2. Once the image was built, run:
```Podman
  podman run -it \
           -v "$(pwd)/logs:/app/logs" \
           -v "$(pwd)/plots:/app/plots" \
           -v "$(pwd)/data:/app/data" \
           schtoppaper:latest

```
- `v` is used to mount 3 of the project folders (`logs`, `plots`, `data`) into a pod, thus allowing all outputs to be saved; please
  note that `data` folder have to be populated with raw data
- the commands starts a `bash` session, from which scripts can be run as described in section above

#### Requirements

To run the scripts, podman machine has to be configured to have at least 24GB of memory, however some
of the scripts (e.g. `2j2_*`), require more than 64GB of memory.

Computations are CPU parallelised, therefore increasing CPU count would increase the speed of computations.

### Docker
*not tested*

1. Build:
```Docker
  docker build -t schtoppaper:latest -f Dockerfile.julia1-8 .
```
2. Once the image was built, run:
```Docker
  docker run -it \
           -v "$(pwd)/logs:/app/logs" \
           -v "$(pwd)/plots:/app/plots" \
           -v "$(pwd)/data:/app/data" \
           schtoppaper:latest
```

- `v` is used to mount 3 of the project folders (`logs`, `plots`, `data`) into a pod, thus allowing all outputs to be saved; please
  note that `data` folder have to be populated with raw data
- the commands starts a `bash` session, from which scripts can be run as described in section above
