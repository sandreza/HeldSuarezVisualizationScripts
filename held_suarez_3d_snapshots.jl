include("utils.jl")
filename = "data/small_earth_snapshots.jld2"
 #filename = "data/earth_snapshots.jld2"
jl_file = jldopen(filename, "r+")

# anomaly
# state_val = jl_file["T"]["1"][:,:,:,1] # ./ jl_file["ρ"]["1"][:,:,:,1]
# state_val = state_val .- mean(state_val, dims = 1)
# velocity
state_val = jl_file["ρw"]["1"][:,:,:,1] ./ jl_file["ρ"]["1"][:,:,:,1]
title_string = "Vertical Velocity"
height_index_start = 1
height_index_end = 31

clims = quantile(state_val[:,:,height_index_start:height_index_end][:], 0.99)
clims = (-clims,clims)

fig = Figure(resolution =(1520, 980))
ax = LScene(fig, scenekw = (camera = cam3d!, show_axis = true))
ax_text = Label(fig, title_string,
        textsize = 30, color = (:black, 0.85))

cmap = :balance # :Blues_9
cmapa = RGBAf0.(to_colormap(cmap), 1);
cmap = vcat(cmapa[1:15], fill(RGBAf0(0,0,0,0), 10), cmapa[25:end])

v1 = volume!(ax, 0..20, 0..10, 0..5, state_val[:,30:180-30,height_index_start:height_index_end],
                 colorrange=clims, algorithm=:absorption, absorption=10.0f0, 
                 colormap=cmap)
axis = ax.scene[OldAxis]
axis[:names, :axisnames] = ("longitude [ᵒ]", "latitude [ᵒ]", "height [km]")
tstyle = axis[:names] #  get the nested attributes and work directly with them

tstyle[:textsize] = 15
tstyle[:textcolor] = (:black, :black, :black)
tstyle[:font] = "helvetica"
tstyle[:gap] = 10
axis[:ticks][:textcolor] = :black
axis[:ticks][:textsize] = 10
cbar1 = Colorbar(fig, v1, label = L" $w$ [m/s]",width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize=25, tickalign = 1, height = Relative(3/4)
)

axis[:ticks][:ranges] = ([0.0, 5.0, 10.0, 15.0, 20.0], [0.0, 2.5, 5.0, 7.5, 10.0], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
axis[:ticks][:labels] = (["180W", "90W", "0", "90E", "180E"], ["60S", "30S", "0", "30N", "60N"], ["0", "6", "12", "18", "24", "30"])

fig[2:10, 1:10] = ax
fig[3:8, 11] = cbar1
fig[1, 5:6] = ax_text

zoom!(ax.scene, 0.8)
update!(fig.scene)

#=
# use this to record video
# bring up the figure and then run this piece of code
seconds = 10
fps = 30
frames = round(Int, fps * seconds )
record(fig, pwd() * "/zonal_velocity.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
=#

##
using LinearAlgebra
# positions
x⃗ = [-1, 0, 1] .* 1.0 

# position operator
x̂ = Diagonal(x⃗)

# construct differentiation matrix / position operator
X⃗ = [x⃗[i]^(j-1) for i in eachindex(x⃗), j in eachindex(x⃗)]
dX⃗ = [(j-1) * x⃗[i]^(j-2) for i in eachindex(x⃗), j in eachindex(x⃗)]
dX⃗[:,1] .= 0.0 # get rid of NaNs
# p̂ * X⃗ = dX⃗, so just invert to find p̂
p̂ = dX⃗ / X⃗

# fake identity
# won't be identiy for polynomial orders ≥ length(x⃗)-1
I_s =  p̂ * x̂ - x̂ * p̂

# and check that its correct to machine precision
check = 3 .+ x⃗ 
I_s * check - check