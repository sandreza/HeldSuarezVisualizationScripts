using GaussQuadrature

function cubed_sphere_warp(
    a,
    b,
    c,
    R = max(abs(a), abs(b), abs(c)),
)

    function f(sR, ξ, η)
        X, Y = tan(π * ξ / 4), tan(π * η / 4)
        ζ1 = sR / sqrt(X^2 + Y^2 + 1)
        ζ2, ζ3 = X * ζ1, Y * ζ1
        ζ1, ζ2, ζ3
    end

    fdim = argmax(abs.((a, b, c)))
    if fdim == 1 && a < 0
        # (-R, *, *) : formulas for Face I from Ronchi, Iacono, Paolucci (1996)
        #              but for us face II of the developed net of the cube
        x1, x2, x3 = f(-R, b / a, c / a)
    elseif fdim == 2 && b < 0
        # ( *,-R, *) : formulas for Face II from Ronchi, Iacono, Paolucci (1996)
        #              but for us face III of the developed net of the cube
        x2, x1, x3 = f(-R, a / b, c / b)
    elseif fdim == 1 && a > 0
        # ( R, *, *) : formulas for Face III from Ronchi, Iacono, Paolucci (1996)
        #              but for us face IV of the developed net of the cube
        x1, x2, x3 = f(R, b / a, c / a)
    elseif fdim == 2 && b > 0
        # ( *, R, *) : formulas for Face IV from Ronchi, Iacono, Paolucci (1996)
        #              but for us face I of the developed net of the cube
        x2, x1, x3 = f(R, a / b, c / b)
    elseif fdim == 3 && c > 0
        # ( *, *, R) : formulas for Face V from Ronchi, Iacono, Paolucci (1996)
        #              and the same for us on the developed net of the cube
        x3, x2, x1 = f(R, b / c, a / c)
    elseif fdim == 3 && c < 0
        # ( *, *,-R) : formulas for Face VI from Ronchi, Iacono, Paolucci (1996)
        #              and the same for us on the developed net of the cube
        x3, x2, x1 = f(-R, b / c, a / c)
    else
        error("invalid case for cubed_sphere_warp(::EquiangularCubedSphere): $a, $b, $c")
    end

    return x1, x2, x3
end

N = 5
ξ¹, ω¹ = GaussQuadrature.legendre(N, GaussQuadrature.both)
ξ² = copy(ξ¹)
ξ³ = copy(ξ¹)

Ne = 12
vertices = collect(range(-1, 1, length = Ne+1))
x¹ = zeros(N, Ne)
x² = copy(x¹)
x³ = copy(x¹)

for e in 1:Ne
    @. x¹[:, e] = (ξ¹ + 1) / 2 * (vertices[e+1] - vertices[e]) + vertices[e]
    @. x²[:, e] = (ξ¹ + 1) / 2 * (vertices[e+1] - vertices[e]) + vertices[e]
    @. x³[:, e] = (ξ¹ + 1) / 2 * (vertices[e+1] - vertices[e]) + vertices[e]
end

face1 = [(ξ¹[1], ξ²[j], ξ³[k]) for j in 1:N, k in 1:N]
warpedface1 = [cubed_sphere_warp(face1[i]...) for i in eachindex(face1)]

face2 = [(ξ¹[end], ξ²[j], ξ³[k]) for j in 1:N, k in 1:N]
warpedface2 = [cubed_sphere_warp(face2[i]...) for i in eachindex(face2)]

face3 = [(ξ¹[j], ξ²[1], ξ³[k]) for j in 1:N, k in 1:N]
warpedface3 = [cubed_sphere_warp(face3[i]...) for i in eachindex(face1)]

face4 = [(ξ¹[j], ξ²[end], ξ³[k]) for j in 1:N, k in 1:N]
warpedface4 = [cubed_sphere_warp(face4[i]...) for i in eachindex(face1)]

face5 = [(ξ¹[j], ξ²[k], ξ³[1]) for j in 1:N, k in 1:N]
warpedface5 = [cubed_sphere_warp(face5[i]...) for i in eachindex(face1)]

face6 = [(ξ¹[j], ξ²[k], ξ³[end]) for j in 1:N, k in 1:N]
warpedface6 = [cubed_sphere_warp(face6[i]...) for i in eachindex(face1)]


##
# this is correct
Nex = Ney = Ne
face1 = [(x¹[1], x²[j, ex], x³[k, ey]) for j in 1:N, k in 1:N, ex in 1:Nex, ey in 1:Ney]
warpedface1 = [cubed_sphere_warp(face1[i]...) for i in eachindex(face1)]

face2 = [(x¹[end], x²[j, ex], x³[k, ey]) for j in 1:N, k in 1:N, ex in 1:Nex, ey in 1:Ney]
warpedface2 = [cubed_sphere_warp(face2[i]...) for i in eachindex(face2)]

face3 = [(x¹[j, ey], x²[1], x³[k, ex]) for j in 1:N, k in 1:N, ex in 1:Nex, ey in 1:Ney]
warpedface3 = [cubed_sphere_warp(face3[i]...) for i in eachindex(face1)]

face4 = [(x¹[j, ey], x²[end], x³[k, ex]) for j in 1:N, k in 1:N, ex in 1:Nex, ey in 1:Ney]
warpedface4 = [cubed_sphere_warp(face4[i]...) for i in eachindex(face1)]

face5 = [(x¹[j, ex], x²[k, ey], x³[1]) for j in 1:N, k in 1:N, ex in 1:Nex, ey in 1:Ney]
warpedface5 = [cubed_sphere_warp(face5[i]...) for i in eachindex(face1)]

face6 = [(x¹[j, ex], x²[k, ey], x³[end]) for j in 1:N, k in 1:N, ex in 1:Nex, ey in 1:Ney]
warpedface6 = [cubed_sphere_warp(face6[i]...) for i in eachindex(face1)]
##
using GLMakie
fullsphere = false
ms = 25
fig, ax, sc = scatter(warpedface1, color = :red, markersize = ms, show_axis = false)
scatter!(ax, warpedface1, color = :red, markersize = ms,)
if fullsphere
    scatter!(ax, warpedface2, color = :red, markersize = ms,)
end
scatter!(ax, warpedface3, color = :blue, markersize = ms,)
if fullsphere
    scatter!(ax, warpedface4, color = :blue, markersize = ms,)
end
# scatter!(ax, warpedface5, color = :green, markersize = ms,)
if fullsphere
    scatter!(ax, warpedface6, color = :green, markersize = ms,)
end
#=
face1 = [(x¹[1], x²[j, ex], x³[k, ey]) for j in 1:N, k in 1:N, ex in 1:Nex, ey in 1:Ney]
warpedface1 = [cubed_sphere_warp(face1[j, k, ex, ey]...) for j in 1:N, k in 1:N, ex in 1:Nex, ey in 1:Ney]

face1_edges1 = [(x¹[1], x²[j, ex], x³[k, ey]) for j in [1, N], k in 1:N, ex in 1:Nex, ey in 1:Ney]
warpedface1_edges1 = [cubed_sphere_warp(face1_edges1[j]...) for j in eachindex(face1_edges1)]

face1_edges2 = [(x¹[1], x²[j, ex], x³[k, ey]) for k in [1, N], j in 1:N, ex in 1:Nex, ey in 1:Ney]
warpedface1_edges2 = [cubed_sphere_warp(face1_edges2[j]...) for j in eachindex(face1_edges2)]

warpedface1 = reshape(warpedface1, (25, 100))

fig, ax, sc = scatter(warpedface1[:, 1], markersize = 0.1, show_axis = false)

for e in 1:Ne^2
    colors = [:red, :green, :orange, :blue, :purple, :yellow, :black, :grey]
    colorindex = (e - 1) % length(colors) + 1 # rand(collect(1:length(colors)))
    colorindex = 1
    scatter!(ax, warpedface1[:, e], color = colors[colorindex], markersize = 20)
end

scatter!(ax, warpedface1_edges1, color = :black, markersize = 20)
scatter!(ax, warpedface1_edges2, color = :black, markersize = 20)
=#


r = range(-1, 1, length = Ne + 1)

xa = zeros(Ne + 1, Ne + 1)
xb = zeros(Ne + 1, Ne + 1)
xc = zeros(Ne + 1, Ne + 1)

xa .= r
xb .= r'
xc .= 1

a = [cubed_sphere_warp(a, b, c)[1] for (a, b, c) in zip(xa, xb, xc)] .* 1.01
b = [cubed_sphere_warp(a, b, c)[2] for (a, b, c) in zip(xa, xb, xc)] .* 1.01
c = [cubed_sphere_warp(a, b, c)[3] for (a, b, c) in zip(xa, xb, xc)] .* 1.01
lw = 1
wireframe!(ax,
    a, b, c,
    show_axis = false,
    linewidth = lw)
#=
wireframe!(sc,
    a, b, -c,
    show_axis = false,
    linewidth = lw)
    =#
#=
wireframe!(sc,
    c, a, b,
    show_axis = false,
    linewidth = lw)
    =#

wireframe!(sc,
    -c, a, b,
    show_axis = false,
    linewidth = lw)
#=
wireframe!(sc,
    b, c, a,
    show_axis = false,
    linewidth = lw)
wireframe!(sc,
    b, -c, a,
    show_axis = false,
    linewidth = lw)
=#

#=
for e in eachindex(warpedface1_edges1)
    scatter!(ax, warpedface1_edges1[e], color = :black, markersize = 20)
    scatter!(ax, warpedface1_edges2[e], color = :black, markersize = 20)
end
=#
#=
r = range(-1,1,length=n+1)

xa = zeros(n+1,n+1)
xb = zeros(n+1,n+1)
xc = zeros(n+1,n+1)

xa .= r
xb .= r'
xc .= 1

a = [cubedshellwarp(a,b,c)[1] for (a,b,c) in zip(xa,xb,xc)]
b = [cubedshellwarp(a,b,c)[2] for (a,b,c) in zip(xa,xb,xc)]
c = [cubedshellwarp(a,b,c)[3] for (a,b,c) in zip(xa,xb,xc)]


lw = 10
wireframe!(sc,
           a,b,c,
           show_axis=false,
           linewidth=lw)
wireframe!(sc,
           a,b,-c,
           show_axis=false,
           linewidth=lw)
wireframe!(sc,
           c,a,b,
           show_axis=false,
           linewidth=lw)
wireframe!(sc,
           -c,a,b,
           show_axis=false,
           linewidth=lw)
wireframe!(sc,
           b,c,a,
           show_axis=false,
           linewidth=lw)
wireframe!(sc,
           b,-c,a,
           show_axis=false,
           linewidth=lw)
=#

##
ms = 25
fig, ax, sc = scatter(warpedface1, color = :red, markersize = 0.0, show_axis = false)

r = range(-1, 1, length = Ne + 1)

xa = zeros(Ne + 1, Ne + 1)
xb = zeros(Ne + 1, Ne + 1)
xc = zeros(Ne + 1, Ne + 1)

xa .= r
xb .= r'
xc .= 1

a = [cubed_sphere_warp(a, b, c)[1] for (a, b, c) in zip(xa, xb, xc)] .* 1.00
b = [cubed_sphere_warp(a, b, c)[2] for (a, b, c) in zip(xa, xb, xc)] .* 1.00
c = [cubed_sphere_warp(a, b, c)[3] for (a, b, c) in zip(xa, xb, xc)] .* 1.00

wireframe!(ax,
    a, b, c,
    show_axis = false,
    linewidth = lw)

wireframe!(sc,
    -c, a, b,
    show_axis = false,
    linewidth = lw*5)

scatter!(ax, warpedface1, color = :red, markersize = ms*0.75,)
scatter!(ax, warpedface3, color = :blue, markersize = ms,)
