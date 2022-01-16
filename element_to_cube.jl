using GaussQuadrature, LinearAlgebra

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

    R  = (ξ³ + 1) / 2 * (r² - r¹) + r¹
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

N = 12
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

vec = ξ¹
ξp = [vec[i]^(j - 1) for i = 1:length(vec), j = 1:length(vec)]
dξp = [(j - 1) * vec[i]^(j - 2) for i = 1:length(vec), j = 1:length(vec)]
∂ξ¹ = dξp / ξp
∂ξ² = copy(∂ξ¹)
∂ξ³ = copy(∂ξ¹)

ω¹ = reshape(ω¹, (N, 1, 1))
ω² = reshape(ω¹, (1, N, 1))
ω³ = reshape(ω¹, (1, 1, N))

##
using GLMakie
fig, ax, sc = scatter(cuby)


##
# Check Metric terms
scale = 1
rectangle = (; aˣ = -1.0/scale, bˣ = 1.0/scale, aʸ = -1.0/scale, bʸ=1.0/scale)
face1 = cuby_sphere.(ξ¹, ξ², ξ³, face = 1, rectangle = rectangle)
L, M, N = size(face1)
x¹ = [face1[i, j, k][1] for i = 1:L, j = 1:M, k = 1:N]
x² = [face1[i, j, k][2] for i = 1:L, j = 1:M, k = 1:N]
x³ = [face1[i, j, k][3] for i = 1:L, j = 1:M, k = 1:N]

∂x¹∂ξ¹ = copy(x¹) .* 0.0
∂x¹∂ξ² = copy(x¹) .* 0.0
∂x¹∂ξ³ = copy(x¹) .* 0.0

∂x²∂ξ¹ = copy(x¹) .* 0.0
∂x²∂ξ² = copy(x¹) .* 0.0
∂x²∂ξ³ = copy(x¹) .* 0.0

∂x³∂ξ¹ = copy(x¹) .* 0.0
∂x³∂ξ² = copy(x¹) .* 0.0
∂x³∂ξ³ = copy(x¹) .* 0.0

jacobian = zeros((3, 3, size(x¹)...))
ijacobian = copy(jacobian)
detjacobian = copy(x¹) .* 0.0

for i = 1:L, j = 1:M, k = 1:N
    for ii = 1:L
        ∂x¹∂ξ¹[i, j, k] += ∂ξ¹[i, ii] * x¹[ii, j, k]
        ∂x²∂ξ¹[i, j, k] += ∂ξ¹[i, ii] * x²[ii, j, k]
        ∂x³∂ξ¹[i, j, k] += ∂ξ¹[i, ii] * x³[ii, j, k]
    end
    for jj = 1:M
        ∂x¹∂ξ²[i, j, k] += ∂ξ²[j, jj] * x¹[i, jj, k]
        ∂x²∂ξ²[i, j, k] += ∂ξ²[j, jj] * x²[i, jj, k]
        ∂x³∂ξ²[i, j, k] += ∂ξ²[j, jj] * x³[i, jj, k]
    end
    for kk = 1:N
        ∂x¹∂ξ³[i, j, k] += ∂ξ³[k, kk] * x¹[i, j, kk]
        ∂x²∂ξ³[i, j, k] += ∂ξ³[k, kk] * x²[i, j, kk]
        ∂x³∂ξ³[i, j, k] += ∂ξ³[k, kk] * x³[i, j, kk]
    end
    jacobian[1, 1, i, j, k] = ∂x¹∂ξ¹[i, j, k]
    jacobian[2, 1, i, j, k] = ∂x²∂ξ¹[i, j, k]
    jacobian[3, 1, i, j, k] = ∂x³∂ξ¹[i, j, k]

    jacobian[1, 2, i, j, k] = ∂x¹∂ξ²[i, j, k]
    jacobian[2, 2, i, j, k] = ∂x²∂ξ²[i, j, k]
    jacobian[3, 2, i, j, k] = ∂x³∂ξ²[i, j, k]

    jacobian[1, 3, i, j, k] = ∂x¹∂ξ³[i, j, k]
    jacobian[2, 3, i, j, k] = ∂x²∂ξ³[i, j, k]
    jacobian[3, 3, i, j, k] = ∂x³∂ξ³[i, j, k]

    detjacobian[i, j, k] = det(jacobian[:, :, i, j, k])
    ijacobian[:,:, i,j,k] = inv(jacobian[:,:, i,j,k])
end

if scale == 1
    volume_numerical = sum(detjacobian .* ω¹ .* ω².* ω³)
    volume_exact = 4/3*π*(1^3 - 0.8^3)/6 # 1/6 since one of the faces
    volume_error = abs(volume_numerical - volume_exact) / abs(volume_exact)
    println("The relative error in computing the volume is ", volume_error)
end

##
# Check Free-Stream
divu = zeros(size(x¹))
c = [1, 1, 1] # vector field in cartesian space
# rows of ijacobian[:,:, i,j,k] are the contravariant basis vectors
for i = 1:L, j = 1:M, k = 1:N
    for ii = 1:L
        u¹ = dot(ijacobian[1, :, ii, j, k], c)
        divu[i, j, k] += ∂ξ¹[i, ii] * detjacobian[ii, j, k] * u¹
    end
    for jj = 1:M
        u² = dot(ijacobian[2, :, i, jj, k], c)
        divu[i, j, k] += ∂ξ²[j, jj] * detjacobian[i, jj, k] * u²
    end
    for kk = 1:N
        u³ = dot(ijacobian[3, :, i, j, kk], c)
        divu[i, j, k] += ∂ξ³[k, kk] * detjacobian[i, j, kk] * u³
    end
end
divu = divu ./ detjacobian
maximum(abs.(divu))


