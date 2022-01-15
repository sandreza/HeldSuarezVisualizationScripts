include("utils.jl")

filename = "data/earth_snapshots.jld2"
# filename = "data/small_earth_snapshots.jld2"
jl_file = jldopen(filename, "r+")

lat_grid = jl_file["grid"]["latitude"]
lon_grid = jl_file["grid"]["longitude"]
rad_grid = jl_file["grid"]["radius"]
λ = lon_grid 
ϕ = lat_grid
r = rad_grid 

fig = Figure(resolution = (1700, 1000))
s_string = "T" # state string
state_file =  jl_file[s_string]["1"][:,:,1,1]
clims = quantile.(Ref(state_file[:]), (0.0,1.0))
colormap =  :afmhot # :thermal # :binary
titlestring =  "Temperature"
ax1 = fig[1,1] = Axis(fig, title = titlestring, titlesize = 30) 
# fig[3,2+3*2] = Label(fig, "ρv: latlon ", textsize = 30) 
hm = heatmap!(ax1, λ, ϕ, state_file, colormap = colormap, colorrange = clims, interpolate = true, shading = false, show_axis=false)
cbar1 = Colorbar(fig, hm, label = L" $T$ [K]",width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize=25, tickalign = 1, height = Relative(3/4)
)
fig[1,2] = cbar1


s_string = "p" # state string
state_file =  jl_file[s_string]["1"][:,:,1,1] ./ 1000
clims = quantile.(Ref(state_file[:]), (0.0,1.0))
colormap = :bone_1 # :acton # :bamako # :balance # :thermal
titlestring =  "Pressure"
ax2 = fig[2,1] = Axis(fig, title = titlestring, titlesize = 30) 
# fig[3,2+3*2] = Label(fig, "ρv: latlon ", textsize = 30) 
hm2 = heatmap!(ax2, λ, ϕ, state_file, colormap = colormap, colorrange = clims, interpolate = true, shading = false, show_axis=false)
cbar2 = Colorbar(fig, hm2, label = L" $P$ [kpa]",width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize=25, tickalign = 1, height = Relative(3/4)
)
fig[2,2] = cbar2

s_string = "ρu" # state string
state_file =  jl_file[s_string]["1"][:,:,1,1] ./ jl_file["ρ"]["1"][:,:,1,1]
clims = quantile.(Ref(state_file[:]), (0.0,1.0))
clims = (-clims[2], clims[2])
colormap = :balance # :thermal
titlestring =  "Zonal Velocity"
ax3 = fig[1,3] = Axis(fig, title = titlestring, titlesize = 30) 
# fig[3,2+3*2] = Label(fig, "ρv: latlon ", textsize = 30) 
hm3 = heatmap!(ax3, λ, ϕ, state_file, colormap = colormap, colorrange = clims, interpolate = true, shading = false, show_axis=false)
cbar3 = Colorbar(fig, hm3, label = L" $u$ [m/s]",width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize=25, tickalign = 1, height = Relative(3/4)
)
fig[1,4] = cbar3

s_string = "ρv" # state string
state_file =  jl_file[s_string]["1"][:,:,1,1] ./ jl_file["ρ"]["1"][:,:,1,1]
clims = quantile.(Ref(state_file[:]), (0.0,1.0))
clims = (-clims[2], clims[2])
colormap = :balance # :thermal
titlestring =  "Meriodonal Velocity"
ax4 = fig[2,3] = Axis(fig, title = titlestring, titlesize = 30) 
# fig[3,2+3*2] = Label(fig, "ρv: latlon ", textsize = 30) 
hm4 = heatmap!(ax4, λ, ϕ, state_file, colormap = colormap, colorrange = clims, interpolate = true, shading = false, show_axis=false)
cbar4 = Colorbar(fig, hm4, label = L" $v$ [m/s]",width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize=25, tickalign = 1, height = Relative(3/4)
)
fig[2,4] = cbar4



for ax in [ax1, ax2, ax3, ax4]
    ax.limits = (extrema(λ)..., extrema(ϕ)...)

    ax.xlabel = "Longitude [ᵒ]"
    ax.ylabel = "Latitude [ᵒ]"
    ax.xlabelsize = 25
    ax.ylabelsize = 25 
    ax.yticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])
    ax.xticks = ([-160, -120, -80, -40, 0, 40, 80, 120, 160], ["160W", "120W", "80W", "40W", "0", "40E", "80E", "120E", "160E"])
end
update!(fig.scene)
display(fig)
