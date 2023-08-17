using JLD2
using Oceananigans.Architectures: arch_array, AbstractArchitecture

import Adapt: adapt_structure, adapt

const non_leep_year_month_days = [15.5, 45.0, 74.5, 105.0, 135.5, 166.0, 196.5, 227.5, 258.0, 288.5, 319.0, 349.5]

struct SurfacePAR{D, I, F, T}
    data :: D

    is :: I
    js :: I

    lon_offset :: F
    lat_offset :: F

    times :: T

    SurfacePAR(data::D, is::I, js::J, lon_offset::F, lat_offset::F, times::T) where {D, I, F, T} = 
        new{D, I, F, T}(data, is, js, lon_offset, lat_offset, times)
end

adapt_structure(to, PAR::SurfacePAR) = SurfacePAR(adapt(to, PAR.data))

@inline function (PAR::SurfacePAR)(x, y, t)
    # bit silly since in the PAR integraiton were going `x = xnode ...`, maybe I should make a discrete version of this
    i, j = unsafe_trunc(Int, x + lon_offset), unsafe_trunc(Int, y + lat_offset)

    i, j = max(1, i), max(1, j)
    i, j = min(is, i), min(js, j)

    return cyclic_interpolate(PAR.data, i, j, t, PAR.times)
end

# required for boundary conditions
@inline function (PAR::SurfacePAR)(i, j, grid, clock, fields)
    t = clock.time

    i, j = max(1, i), max(1, j)
    i, j = min(is, i), min(js, j)

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
                    times = non_leep_year_month_days)

    surfac_PAR_file = jldopen(data_path)

    surfac_PAR_data = surfac_PAR_file[variable_name]
    surfac_PAR_data[isnan.(surfac_PAR_data)] .= 0.0

    surfac_PAR_data = arch_array(architecture, surfac_PAR_data)

    close(surfac_PAR_file)

    return SurfacePAR(surfac_PAR_data, is, js, lon_offset, lat_offset, times)
end

