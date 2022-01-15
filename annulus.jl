using GaussQuadrature, LinearAlgebra

# First construct the mapping for a spherical sector

# parameters:
nn = 12  # number of sectors
φ = π / nn # half angle of sector
r1 = 0.8 # inner radius
r2 = 1.0 # outer radius

# mappings
# annulus mappings
R(ξ¹; r1 = r1, r2 = r2) = (r2 - r1) * (ξ¹ + 1) / 2 + r1
Θ(ξ²; θ1 = π / 2 - φ, θ2 = π / 2 + φ) = (θ2 - θ1) * (ξ² + 1) / 2 + θ1

# annulus to cartesian mappings
x_mapping(r, θ) = r * cos(θ)
y_mapping(r, θ) = r * sin(θ)

# Now construct the discretization
# Construct Gauss-Lobatto points
N1 = 1
ξ1vec, ω1 = GaussQuadrature.legendre(N1 + 1, GaussQuadrature.both)
ω1 = reshape(ω1, (N1 + 1, 1))
N2 = 2
ξ2vec, ω2 = GaussQuadrature.legendre(N2 + 1, GaussQuadrature.both)
ω2 = reshape(ω2, (1, N2 + 1))

# Construct Differentiation Matrices
vec = ξ1vec
ξp = [vec[i]^(j - 1) for i in 1:length(vec), j in 1:length(vec)]
dξp = [(j - 1) * vec[i]^(j - 2) for i in 1:length(vec), j in 1:length(vec)]
∂ξ¹ = dξp / ξp

vec = ξ2vec
ξp = [vec[i]^(j - 1) for i in 1:length(vec), j in 1:length(vec)]
dξp = [(j - 1) * vec[i]^(j - 2) for i in 1:length(vec), j in 1:length(vec)]
∂ξ² = dξp / ξp

# Calculate the grid points in polar coordinates
rvec = R.(ξ1vec)
θvec = Θ.(ξ2vec)
# convert to cartesian
x_positions = [x_mapping(r, θ) for r in rvec, θ in θvec]
y_positions = [y_mapping(r, θ) for r in rvec, θ in θvec]

# Now do the discrete Calculus
∂x∂ξ¹ = copy(x_positions)
∂y∂ξ¹ = copy(x_positions)
for j in 1:size(∂x∂ξ¹)[2]
    ∂x∂ξ¹[:, j] = ∂ξ¹ * x_positions[:, j]
    ∂y∂ξ¹[:, j] = ∂ξ¹ * y_positions[:, j]
end

∂x∂ξ² = copy(x_positions)
∂y∂ξ² = copy(x_positions)
for i in 1:size(∂x∂ξ²)[1]
    ∂x∂ξ²[i, :] = ∂ξ² * x_positions[i, :]
    ∂y∂ξ²[i, :] = ∂ξ² * y_positions[i, :]
end
jacobian = zeros(size(x_positions)..., (2, 2)...)

# Construct convenient object for entries, The columns are the covariant vectors
jacobian[:, :, 1, 1] .= ∂x∂ξ¹
jacobian[:, :, 1, 2] .= ∂x∂ξ²
jacobian[:, :, 2, 1] .= ∂y∂ξ¹
jacobian[:, :, 2, 2] .= ∂y∂ξ²

detJ = [det(jacobian[i, j, :, :]) for i in 1:length(ξ1vec), j in 1:length(ξ2vec)]
M = detJ .* ω1 .* ω2
exact_area = (r2^2 - r1^2) * φ
approx_area = sum(M)

wrongness = (exact_area - approx_area) / exact_area
println("The error in computing the area is ", wrongness)

# Construct contravariant basis numerically, these are the face normals
ijacobian = copy(jacobian) # the rows are the contravariant vectors
for i in 1:length(ξ1vec), j in 1:length(ξ2vec)
    tmp = inv(jacobian[i, j, :, :])
    ijacobian[i, j, :, :] .= tmp
end
ijacobian[1, 1, :, :]
x_positions[1, 1], y_positions[1, 1]

# face 2 is the linear side
approx_vec = ijacobian[1, 1, 2, :] ./ norm(ijacobian[1, 1, 2, :])
exact_vec = [-y_positions[1, 1], x_positions[1, 1]] ./ norm([y_positions[1, 1], -x_positions[1, 1]])
println("angle face error ", norm(approx_vec - exact_vec) / norm(exact_vec))
# face 1 is the curvy side
approx_vec = ijacobian[1, 1, 1, :] ./ norm(ijacobian[1, 1, 1, :])
exact_vec = [x_positions[1, 1], y_positions[1, 1]] ./ norm([x_positions[1, 1], y_positions[1, 1]])
println("radial face error ", norm(approx_vec - exact_vec) / norm(exact_vec))

# Analytic Metrics 
M, N = size(detJ)
a₁ᴬ = zeros(M, N, 2)
a₂ᴬ = zeros(M, N, 2)
jacobianᴬ = zeros(M, N, 2, 2)
ijacobianᴬ = zeros(M, N, 2, 2)
detJᴬ = zeros(M, N)
a¹ᴬ = zeros(M, N, 2)
a²ᴬ = zeros(M, N, 2)
for i in 1:M, j in 1:N
    a₁ᴬ[i, j, 1] = (rvec[2] - rvec[1]) / 2 * cos((θvec[end] - θvec[1]) * (ξ2vec[j] + 1) / 2 + θvec[1])
    a₁ᴬ[i, j, 2] = (rvec[2] - rvec[1]) / 2 * sin((θvec[end] - θvec[1]) * (ξ2vec[j] + 1) / 2 + θvec[1])

    a₂ᴬ[i, j, 1] = -((rvec[2] - rvec[1]) / 2 * (ξ1vec[i] + 1) + rvec[1]) * sin((θvec[end] - θvec[1]) * (ξ2vec[j] + 1) / 2 + θvec[1]) * (θvec[end] - θvec[1]) / 2
    a₂ᴬ[i, j, 2] = ((rvec[2] - rvec[1]) / 2 * (ξ1vec[i] + 1) + rvec[1]) * cos((θvec[end] - θvec[1]) * (ξ2vec[j] + 1) / 2 + θvec[1]) * (θvec[end] - θvec[1]) / 2

    jacobianᴬ[i, j, 1, 1] = a₁ᴬ[i, j, 1]
    jacobianᴬ[i, j, 2, 1] = a₁ᴬ[i, j, 2]
    jacobianᴬ[i, j, 1, 2] = a₂ᴬ[i, j, 1]
    jacobianᴬ[i, j, 2, 2] = a₂ᴬ[i, j, 2]
    detJᴬ[i, j] = det(jacobianᴬ[i, j, :, :])
    ijacobianᴬ[i, j, :, :] = inv(jacobianᴬ[i, j, :, :])

    a¹ᴬ[i, j, 1] = ijacobianᴬ[i, j, 1, 1]
    a¹ᴬ[i, j, 2] = ijacobianᴬ[i, j, 1, 2]

    a²ᴬ[i, j, 1] = ijacobianᴬ[i, j, 2, 1]
    a²ᴬ[i, j, 2] = ijacobianᴬ[i, j, 2, 2]
end

# Free-Stream Property
i = 1
j = 1
ijacobian[i, j, :, :]
jacobian[i, j, :, :]

a₁ = jacobian[:, :, :, 1];
a₂ = jacobian[:, :, :, 2];

a¹ = ijacobian[:, :, 1, :];
a² = ijacobian[:, :, 2, :];

dot(a¹[1, 1, :], a₁[1, 1, :])
dot(a¹[1, 1, :], a₂[1, 1, :])
dot(a²[1, 1, :], a₁[1, 1, :])
dot(a²[1, 1, :], a₂[1, 1, :])
J = reshape(detJ, (size(a¹)[1:2]..., 1))

# choose vector field in cartesian coordinates
c = [1 1]
# project onto contravariant components
u¹ = [(c*a¹[i, j, :])[1] for i = 1:M, j = 1:N]
u² = [(c*a²[i, j, :])[1] for i = 1:M, j = 1:N]

divu = zeros(size(u¹))
detiJ = 1 ./ detJ
for i = 1:M, j = 1:N
    for ii = 1:M
        divu[i, j] += detiJ[i, j] * ∂ξ¹[i, ii] * detJ[ii, j] * u¹[ii, j]
    end
end
for j = 1:N, i = 1:M
    for jj = 1:N
        divu[i, j] += detiJ[i, j] * ∂ξ²[j, jj] * detJ[i, jj] * u²[i, jj]
    end
end
println("numerically computed divergence with numerical metrics")
display(divu)

# choose vector field in cartesian coordinates
c = [1 1]
# project onto contravariant components
u¹ = [(c*a¹ᴬ[i, j, :])[1] for i = 1:M, j = 1:N]
u² = [(c*a²ᴬ[i, j, :])[1] for i = 1:M, j = 1:N]

divu = zeros(size(u¹))
detiJᴬ = 1 ./ detJᴬ
for i = 1:M, j = 1:N
    for ii = 1:M
        divu[i, j] += detiJᴬ[i, j] * ∂ξ¹[i, ii] * detJᴬ[ii, j] * u¹[ii, j]
    end
end
for j = 1:N, i = 1:M
    for jj = 1:N
        divu[i, j] += detiJᴬ[i, j] * ∂ξ²[j, jj] * detJᴬ[i, jj] * u²[i, jj]
    end
end
println("numerically computed divergence with analytic metrics")
display(divu)










