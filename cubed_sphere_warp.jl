using GaussQuadrature, GLMakie

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
vertices = collect(range(-1, 1, length = Ne + 1))
x¹ = zeros(N, Ne)
x² = copy(x¹)
x³ = copy(x¹)

for e in 1:Ne
    @. x¹[:, e] = (ξ¹ + 1) / 2 * (vertices[e+1] - vertices[e]) + vertices[e]
    @. x²[:, e] = (ξ¹ + 1) / 2 * (vertices[e+1] - vertices[e]) + vertices[e]
    @. x³[:, e] = (ξ¹ + 1) / 2 * (vertices[e+1] - vertices[e]) + vertices[e]
end

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

# Now Plot
lw = 1  # linewidth
ms = 25 # marker size

fig = Figure(resolution = (1400, 1200))
ax = LScene(fig[1, 1])
# fig, ax, sc = scatter(warpedface1, color = :red, markersize = 0.0, show_axis = false)
# scatter(warpedface1, color = :red, markersize = 0.0, show_axis = false)
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

wireframe!(ax,
    -c, a, b,
    show_axis = false,
    linewidth = lw * 5)

scatter!(ax, warpedface1, color = :red, markersize = ms * 0.75,)
scatter!(ax, warpedface3, color = :blue, markersize = ms,)

rotate_cam!(fig.scene.children[1], (π/12, π, 0))
