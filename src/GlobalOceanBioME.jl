module GlobalOceanBioME

export OneDegreeSurfacePAR

using DataDeps, Oceananigans.Units

function __init__(; remove_existing_data=false)

    ## Data for the one_degree_global_simulation
    branch_url = "https://github.com/OceanBioME/OceanBioMEArtifacts/raw/main"
    dir = "near_global/one_degree/2010"
    par_name = "PAR.jld2"
    initial_conditions_name = "initial_conditions.jld2"

    par_url = joinpath(branch_url, dir, par_name)
    initial_conditions_url = joinpath(branch_url, dir, initial_conditions_name)

    dep = DataDep("2010_near_global_bgc",
                  "PAR and intiial conditions for proof of concept near global one degree bgc model ",
                  [par_url, initial_conditions_url])

    DataDeps.register(dep)
    remove_existing_data && rm(datadep"near_global_one_degree_bgc", recursive=true, force=true)
end

include("cyclic_interpolation.jl")
include("surface_par.jl")
include("one_degree_surface_par.jl")

end # module GlobalOceanBioME
