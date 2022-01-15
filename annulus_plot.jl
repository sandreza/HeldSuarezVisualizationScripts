include("annulus.jl")
using GLMakie

lw = 3
Nθ = 100
θ = collect(range(0, 2π, length = Nθ))
xp = cos.(θ)
yp = sin.(θ)
fig = Figure(resolution = (1450, 1350))
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", ylabelsize = 32,
    xlabelsize = 32, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
    xticksize = 30, ytickalign = 1, yticksize = 30, xlabelpadding = -10,
    xticklabelsize = 30, yticklabelsize = 30)
lines!(ax, xp .* r2, yp .* r2, color = :black, linewidth = lw)
lines!(ax, xp .* r1, yp .* r1, color = :black, linewidth = lw)

θs = collect(0:nn) * 2π / nn
for θ in θs
    ray_x = [0.8, 1] * cos(π / 2 + θ + φ)
    ray_y = [0.8, 1] * sin(π / 2 + θ + φ)
    lines!(ax, ray_x, ray_y, color = :black, linewidth = lw)
    ray_x = [0.8, 1] * cos(π / 2 + θ - φ)
    ray_y = [0.8, 1] * sin(π / 2 + θ - φ)
    lines!(ax, ray_x, ray_y, color = :black, linewidth = lw)
end

scatter!(ax, x_positions[2:end], y_positions[2:end], color = :red, markersize = 30)
scatter!(ax, x_positions[1:1], y_positions[1:1], color = :red, markersize = 60, marker = :star5)