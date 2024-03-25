# GlobalOceanBioME
This repository contains code for modelling global ocean biogeochemistry with [OceanBioME](https://github.com/oceanBioME/oceanBioME.jl) and [Oceananigans](https://github.com/cliMA/oceananigans.jl/).

Currently, it contains a script for a repeat year forcing/free runing near global 1 degree model using [ClimaOcean](https://github.com/cliMA/ClimaOcean.jl) as a base. 

**Please note** this is a work in progress and the scripts were written as demonstrations - the physics and biogeochemistry are not complete or calibrated.

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
