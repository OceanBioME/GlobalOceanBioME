using Oceananigans, CairoMakie, GeoMakie

P = FieldTimeSeries("bgc/near_global_360_150_48_ConvectiveAdjustmentVerticalDiffusivity_fields.jld2", "P")

n = Observable(1)

P_plt = @lift interior(P[$n], :, :, 1)

fig = Figure()

ax = GeoAxis(fig[1, 1];title, dest = "+proj=natearth2")

x, y, z = nodes(P)

hm = heatmap!(ax, x, y, P_plt, colorrange = (0.0001, 5), colorscale = log10)

Colorbar(fig[2, 1], hm, vertical= false, label = "Phytoplankton (mmolN/m³)")

CairoMakie.record(fig, "P.mp4", 1:10:length(P.times)) do i
    n[] = i
    @info "$(n.val)"
end
