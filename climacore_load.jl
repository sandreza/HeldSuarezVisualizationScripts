using NCDatasets, GLMakie
using Statistics: mean

datapath = "CHANGE_ME!!!!!"
datafiles = ["held_suarez_rhoe_day" * string(t) * ".nc" for t in 200:10:1190];

datafile = datafiles[1]
ncdata = NCDataset(datapath * datafile, "r")
lon = ncdata["lon"][:]
lat = ncdata["lat"][:]
z = ncdata["z"][:]
time = ncdata["time"][:]
rho = ncdata["rho"][:]
e_tot = ncdata["e_tot"][:]
u = ncdata["u"][:]
v = ncdata["v"][:]

for datafile in datafiles[2:end]
    local ncdata = NCDataset(datapath * datafile, "r")
    global time = append!(time, ncdata["time"][:])
    global rho = cat(rho, ncdata["rho"][:], dims=4)
    global e_tot = cat(e_tot, ncdata["e_tot"][:], dims=4)
    global u = cat(u, ncdata["u"][:], dims=4)
    global v = cat(v, ncdata["v"][:], dims=4)
end

Φ = similar(rho);
for iz in 1:length(z)
    Φ[:, :, iz, :] .= 9.8 * z[iz]
end
K = 0.5 * (u .* u + v .* v);
e_int = e_tot .- K .- Φ;

cv_d = 717.6450000000001;
R_d = 287.058;
cp_d = 1004.7030000000002;
T_0 = 273.16;
MSLP = 1e5;
T = e_int / cv_d .+ T_0;
p = rho * R_d .* T;
θ = T .* (MSLP ./ p) .^ (R_d / cp_d);

#= 
fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, u[:, :, 1, 1], colormap=:balance, interpolate=true, colorrange=(-14, 14))
=#

x = zeros(length(lon), length(lat))
y = zeros(length(lon), length(lat))
z = zeros(length(lon), length(lat))

for (i, loni) in enumerate(lon), (j, latj) in enumerate(lat)
    λ_azimuthal = loni + 180 # λ_azimuthal ∈ [0°, 360°]
    φ_azimuthal = 90 .- latj  # φ_azimuthal ∈ [0°, 180°] (0° at north pole)
    x[i, j] = @. cosd(λ_azimuthal) * sind(φ_azimuthal)
    y[i, j] = @. sind(λ_azimuthal) * sind(φ_azimuthal)
    z[i, j] = @. cosd(φ_azimuthal)
end

itime = Observable(1)
t_iter = 1:length(time)
fig = Figure()
ax = LScene(fig[1, 1])
surface!(ax, x, y, z, color=@lift(T[:, :, 1, $itime]), colormap=:afmhot, colorrange=(256.77658780587706, 299.868732819969))

rotation = (0 * π / 5, π / 6, 0)
rotate_cam!(ax.scene, rotation)

record(fig, "T.mp4", t_iter; framerate=5) do it
    itime[] = it
end