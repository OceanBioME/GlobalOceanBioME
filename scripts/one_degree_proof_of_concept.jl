using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation, cyclic_interpolate, current_time_index, next_time_index

using Oceananigans, OceanBioME, GlobalOceanBioME, DataDeps, JLD2
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    MixingLength, TurbulentKineticEnergyEquation, CATKEVerticalDiffusivity
using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity

using Oceananigans.Architectures: arch_array

import Oceananigans.TurbulenceClosures: cell_diffusion_timescale

import Adapt: adapt_structure, adapt

architecture = CPU()

cell_diffusion_timescale(::IsopycnalSkewSymmetricDiffusivity, args...) = Inf

#####
##### Boundary layer turbulence closure options
#####

convective_κz = 0.1
background_κz = 0.0
background_νz = 0.01
convective_νz = 0.0

cavd = ConvectiveAdjustmentVerticalDiffusivity(; convective_κz, background_κz, background_νz, convective_νz)

boundary_layer_turbulence_closure = cavd

@show boundary_layer_turbulence_closure

#####
##### Build the simulation
#####

start_time = 345days
stop_time = start_time + 2*365days
with_isopycnal_skew_symmetric_diffusivity = true

biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus

simulation = one_degree_near_global_simulation(architecture; 
    start_time, stop_time,
    with_isopycnal_skew_symmetric_diffusivity,
    boundary_layer_turbulence_closure,
    isopycnal_κ_skew = 900.0,
    isopycnal_κ_symmetric = 900.0,
    interior_background_vertical_viscosity = 1e-4,
    surface_background_vertical_viscosity = 1e-4,
    biogeochemistry,
    biogeochemistry_kwargs = (surface_photosynthetically_active_radiation = OneDegreeSurfacePAR(architecture), ),
    tracers = (:N, :P, :Z, :D, :T, :S) # have to specify since NPZD adds T but not S, and buoyancy requires both
)

# Define output
slices_save_interval = 2day
fields_save_interval = 1days
Nx, Ny, Nz = size(simulation.model.grid)

dir = "bgc" 
closure_name = typeof(boundary_layer_turbulence_closure).name.wrapper
output_prefix = "near_global_$(Nx)_$(Ny)_$(Nz)_$closure_name"
with_isopycnal_skew_symmetric_diffusivity || (output_prefix *= "_no_gm")

simulation.output_writers[:checkpointer] = Checkpointer(simulation.model; dir,
                                                        prefix = output_prefix * "_checkpointer",
                                                        schedule = TimeInterval(365days),
                                                        cleanup = true,
                                                        overwrite_existing = true)

model = simulation.model

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers, (PAR = model.biogeochemistry.light_attenuation_model.field, )); dir,
                                                      schedule = TimeInterval(fields_save_interval),
                                                      filename = output_prefix * "_fields",
                                                      with_halos = true,
                                                      overwrite_existing = true)

""" Load initial conditions from Copernicus models.

P and N are a direct downsampling and unit conversion from https://doi.org/10.48670/moi-00015
Z is downsampled and unit converted from https://doi.org/10.48670/moi-00020 and then divided by 
the mixed layer depth from https://doi.org/10.48670/moi-00024, and then applied uniformly over
the mixed region
"""
file = jldopen(datadep"2010_near_global_bgc/initial_conditions.jld2")

N_init = arch_array(architecture, file["N"])
P_init = arch_array(architecture, file["P"])
Z_init = arch_array(architecture, file["Z"])

close(file)

set!(model, N = N_init, P = P_init, Z = Z_init)

scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:N, :P, :Z, :D))
simulation.callbacks[:nan_tendencies] = Callback(remove_NaN_tendencies!; callsite = TendencyCallsite())

simulation.callbacks[:nan_checker] = Callback(Oceananigans.Simulations.NaNChecker(; fields = merge(model.tracers, model.velocities), erroring = true), IterationInterval(1))

simulation.Δt = 2minute
simulation.stop_time = time(simulation) + 1days

@info "Running a simulation with Δt = $(prettytime(simulation.Δt)) from $(prettytime(simulation.model.clock.time)) until $(prettytime(simulation.stop_time))"

run!(simulation)

simulation.callbacks[:nan_checker] = Callback(Oceananigans.Simulations.NaNChecker(; fields = merge(model.tracers, model.velocities), erroring = true), IterationInterval(10))

simulation.Δt = 20minutes
simulation.stop_time = start_time + 3 * 365days

@info "Running a simulation with Δt = $(prettytime(simulation.Δt)) from $(prettytime(simulation.model.clock.time)) until $(prettytime(simulation.stop_time))"

run!(simulation; pickup = false)