using GaussQuadrature

function cuby_sphere(ξ¹, ξ², ξ³; radii = nothing, rectangle = nothing, face = 1)
    if isnothing(radii)
        r², r¹ = 1.0, 0.8
    else
        r², r¹ = radii.r², radii.r¹
    end
    if isnothing(rectangle)
        bˣ, aˣ = 1.0, -1.0
        bʸ, aʸ = 1.0, -1.0
    else
        aˣ, bˣ, aʸ, bʸ = rectangle.aˣ, rectangle.bˣ, rectangle.aʸ, rectangle.bʸ
    end

    R = (ξ³ + 1) / 2 * (r² - r¹) + r¹
    η² = (ξ¹ + 1) / 2 * (bˣ - aˣ) + aˣ
    η³ = (ξ² + 1) / 2 * (bʸ - aʸ) + aʸ
    ρ = sqrt(1 + tan(π / 4 * η²)^2 + tan(π / 4 * η³)^2)
    x¹ = R / ρ
    x² = R * tan(π / 4 * η²) / ρ
    x³ = R * tan(π / 4 * η³) / ρ
    if face == 1
        return x¹, x², x³
    elseif face == 2
        return -x¹, -x², -x³
    elseif face == 3
        return x², x¹, x³
    elseif face == 4
        return -x², -x¹, -x³
    elseif face == 5
        return x³, x², x¹
    elseif face == 6
        return -x³, -x², -x¹
    else
        return nothing
    end
    return nothing
end

N = 22
ξ¹, ω¹ = GaussQuadrature.legendre(N, GaussQuadrature.both)
ξ² = copy(ξ¹)
ξ³ = copy(ξ¹)

ξ¹ = reshape(ξ¹, (N, 1, 1))
ξ² = reshape(ξ¹, (1, N, 1))
ξ³ = reshape(ξ¹, (1, 1, N))

cuby = cuby_sphere.(ξ¹, ξ², ξ³)[:]

cuby1 = cuby_sphere.(ξ¹, ξ², ξ³, face = 1)[:]
cuby2 = cuby_sphere.(ξ¹, ξ², ξ³, face = 2)[:]

cuby3 = cuby_sphere.(ξ¹, ξ², ξ³, face = 3)[:]
cuby4 = cuby_sphere.(ξ¹, ξ², ξ³, face = 4)[:]

cuby5 = cuby_sphere.(ξ¹, ξ², ξ³, face = 5)[:]
cuby6 = cuby_sphere.(ξ¹, ξ², ξ³, face = 6)[:]

cuby = [cuby1..., cuby2...]

cuby = [cuby..., cuby3..., cuby4...]
cuby = [cuby..., cuby5..., cuby6...]

##
using GLMakie
scatter(cuby)