using JLD2
using Oceananigans.Architectures: on_architecture, AbstractArchitecture

import Adapt: adapt_structure, adapt

const non_leep_year_month_days = [15.5, 45.0, 74.5, 105.0, 135.5, 166.0, 196.5, 227.5, 258.0, 288.5, 319.0, 349.5]

struct SurfacePAR{D, I, F, T}
    data :: D

    is :: I
    js :: I

    lon_offset :: F
    lat_offset :: F

    times :: T

    resolution :: F

    SurfacePAR(data::D, is::I, js::I, lon_offset::F, lat_offset::F, times::T, resolution::F) where {D, I, F, T} = 
        new{D, I, F, T}(data, is, js, lon_offset, lat_offset, times, resolution)
end

adapt_structure(to, PAR::SurfacePAR) = SurfacePAR(adapt(to, PAR.data),
                                                  adapt(to, PAR.is),
                                                  adapt(to, PAR.js),
                                                  adapt(to, PAR.lon_offset),
                                                  adapt(to, PAR.lat_offset),
                                                  adapt(to, PAR.times),
                                                  adapt(to, PAR.resolution))
@inline function (PAR::SurfacePAR)(x, y, t)
    # bit silly since in the PAR integraiton were going `x = xnode ...`, maybe I should make a discrete version of this
    i, j = unsafe_trunc(Int, (x + PAR.lon_offset) * PAR.resolution), unsafe_trunc(Int, (y + PAR.lat_offset) * PAR.resolution)

    i, j = max(1, i), max(1, j)
    i, j = min(PAR.is, i), min(PAR.js, j)

    return cyclic_interpolate(PAR.data, i, j, t, PAR.times)
end

# required for boundary conditions
@inline function (PAR::SurfacePAR)(i, j, grid, clock, fields)
    t = clock.time

    i, j = max(1, i), max(1, j)
    i, j = min(PAR.is, i), min(PAR.js, j)

    return cyclic_interpolate(PAR.data, i, j, t, PAR.times)
end

# Data origionally from NASA MODIS-Aqua 2010 monthly mean (https://oceandata.sci.gsfc.nasa.gov)
function SurfacePAR(architecture::AbstractArchitecture; 
                    lon_offset = 180.5,
                    lat_offset = 75.5,
                    is = 360,
                    js = 150,
                    data_path = datadep"2010_near_global_bgc/PAR.jld2",
                    variable_name = "one_degree_climatology",
                    times = non_leep_year_month_days,
                    resolution = 1)

    surfac_PAR_file = jldopen(data_path)

    surfac_PAR_data = surfac_PAR_file[variable_name]
    surfac_PAR_data[isnan.(surfac_PAR_data)] .= 0.0

    surfac_PAR_data = on_architecture(architecture, surfac_PAR_data)

    close(surfac_PAR_file)

    return SurfacePAR(surfac_PAR_data, is, js, lon_offset, lat_offset, times, resolution)
end

