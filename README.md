# GlobalOceanBioME
This repository contains code for modelling global ocean biogeochemistry with [OceanBioME](https://github.com/oceanBioME/oceanBioME.jl) and [Oceananigans](https://github.com/cliMA/oceananigans.jl/).

Currently, it contains a script for a repeat year forcing/free runing near global 1 degree model using [ClimaOcean](https://github.com/cliMA/ClimaOcean.jl) as a base. 

**Please note** this is a work in progress and the scripts were written as demonstrations - the physics and biogeochemistry are not complete or calibrated.

https://github.com/OceanBioME/GlobalOceanBioME/assets/26657828/ca2dd956-5fa5-4bb3-b564-743cbd379bd7

<details>
  <summary>More detailed description of proof of concept above</summary>
  This proof of concept is a 1°x1° near global run between 75°N and 75°S with 48 vertical levels.
  
  The physics includes isopycnal skew-symmetric diffusivity and a convective adjustment boundary layer closure, and the biogeochemistry is a simple 4-variable nutrient, phytoplankton, zooplankton, detritus model (<a href="https://oceanbiome.github.io/OceanBioME.jl/stable/model_components/biogeochemical/NPZ/">detailed here</a>). 
  
  The physics is initialised with temperature and salinity, and the top boundary is forced with data, from ECCO v4 (as per <a href="https://github.com/CliMA/ClimaOcean.jl">ClimaOcean NearGlobalSimulation</a>), but the biogeochemistry in this simple example is initialised to be uniform across the globe (hence the abrupt change at the start). The surface photosynthetically available radiation for the biogeochemical model is interpolated from NASA ocean colour data.

  This model runs in around 45 minutes per year to run on an Nvidia A100 GPU.
</details>

## To run the proof of concept script

```
$ git clone git@github.com:OceanBioME/GlobalOceanBioME.git

$ cd GlobalOceanBioME

$ julia --project

julia> ]

(GlobalOceanBioME) > instantiate

julia> include("scripts/one_degree_proof_of_concept.jl")
```

We very much welcome contributors, please contact us if you are interested in collaborating on this project.
