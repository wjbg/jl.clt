using Printf
using DocStringExtensions
using JSON

# Reuter matrix used to rotate stiffness matrix.
global const R = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 2.0]

"""
    Object to represent a material.

### Fields

- name::String
    Identifier or name for the material.
- E1::Float64
    Modulus in 1-direction.
- E2::Float64
    Modulus in 2-direction.
- G12::Float64
    Shear modulus.
- ν12::Float64
     Poisson's ratio (12).
- α1::Float64
    Coefficient of thermal expansion in 1-direction.
- α2::Float64
    Coefficient of thermal expansion in 2-direction.
- S1t::Float64
    Tensile strength in 1-direction.
- S1c::Float64
    Comparessive strength in 1-direction.
- S2t::Float64
    Tensile strength in 2-direction.
- S2c::Float64
    Comparessive strength in 2-direction.
- S6::Float64
    Shear strength.
"""
struct Material
    name::String
    E1::Float64
    E2::Float64
    G12::Float64
    ν12::Float64
    α1::Float64
    α2::Float64
    S1t::Float64
    S1c::Float64
    S2t::Float64
    S2c::Float64
    S6::Float64
end

"""
    Material(d::Dict)

Creates a Material object from a dictionary.

### Arguments

- d::Dict
    Dictionary with the the following fields:
      - name
          Identifier for material.
      - E1, E2, G12, v12
          Moduli and Poisson ratio.
      - alpha1, alpha2
          Coefficients of thermal expansion.
      - S1t, S1c, S2t, S2c, S6
          Strength values.
"""
function Material(d::Dict)
    return Material(d["name"],
                    d["E1"], d["E2"], d["G12"], d["v12"],
                    d["alpha1"], d["alpha2"],
                    d["S1t"], d["S1c"], d["S2t"], d["S2c"], d["S6"])
end

"""
    Material(fn::String)

Creates a Material object from a JSON file.

### Arguments

- fn:String
    Filename of JSON file to load.

The JSON file should have the following fields.
- name
    Identifier for material.
- E1, E2, G12, v12
    Moduli and Poisson ratio.
- alpha1, alpha2
    Coefficients of thermal expansion.
- S1t, S1c, S2t, S2c, S6
    Strength values.
"""
function Material(fn::String)
    d = JSON.parsefile(fn)
    return Material(d)
end

function Base.show(io::IO, mat::Material)
    @printf(io, "%s ", mat.name)
    @printf(io, "(E1: %.2f GPa, ", mat.E1/1E9)
    @printf(io, "E2: %.2f GPa, ", mat.E2/1E9)
    @printf(io, "G12: %.2f GPa, ", mat.G12/1E9)
    @printf(io, "ν12: %.2f, ", mat.ν12)
    @printf(io, "α1: %.2f μm/m, ", mat.α1*1E6)
    @printf(io, "α2: %.2f μm/m)", mat.α2*1E6)
end

"""
Type that represents a composite ply.

### Fields

- mat::Material
    Ply material.
- θ::Float64
    Ply orientation [rad].
- t::Float64
    Ply thickness.
"""
struct Ply
    mat::Material
    θ::Float64
    t::Float64
end

function Base.show(io::IO, ply::Ply)
    @printf(io, "Ply (%s, %.1f°, %.2f mm)", ply.mat.name, rad2deg(ply.θ), ply.t*1e3)
end

"""
    stiffness_matrix(mat::Material)

Compute the stiffness matrix of `mat`.
"""
function stiffness_matrix(mat::Material)
    E1, E2, G12 = mat.E1, mat.E2, mat.G12
    ν12, ν21 = mat.ν12, mat.E2*mat.ν12/mat.E1
    C = [E1/(1-ν12*ν21)     ν21*E1/(1-ν12*ν21) 0;
         ν21*E1/(1-ν12*ν21) E2/(1-ν12*ν21)     0;
         0                  0                  G12]
end

"""
    stiffness_matrix(mat::Material, θ::Float64)

Compute the stiffness matrix of `mat` when when the material coordinate is
rotated by angle `θ` [rad].
"""
function stiffness_matrix(mat::Material, θ::Float64)
    C = stiffness_matrix(mat)
    T = transformation_matrix(θ)
    return inv(T)*C*R*T*inv(R)
end

"""
    stiffness_matrix(ply::Ply)

Compute the stiffness matrix of `ply`.
"""
function stiffness_matrix(ply::Ply)
    return stiffness_matrix(ply.mat, ply.θ)
end

"""
    cte_vector(mat::Material)

Compute the vector with coefficients of thermal expansions of `mat`.
"""
function cte_vector(mat::Material)
    return [mat.α1, mat.α2, 0]
end

"""
    cte_vector(mat::Material, θ::Float64)

Compute the vector with coefficients of thermal expansions of `mat` when the
the coordinate system is rotated by angle `θ` [rad].
"""
function cte_vector(mat::Material, θ::Float64)
    α = cte_vector(mat)
    T = transformation_matrix(θ)
    return R*inv(T)*inv(R)*α
end

"""
    cte_vector(ply::Ply)

Compute the vector with coefficients of thermal expansions of `ply`.
"""
function cte_vector(ply::Ply)
    return cte_vector(ply.mat, ply.θ)
end

"""
    transformation_matrix(θ::Float64)

Compute the transformation matrix for angle `θ` [rad].
"""
function transformation_matrix(θ::Float64)
    n, m = sin(θ), cos(θ)
    T = [m^2  n^2 2*n*m;
         n^2  m^2 -2*n*m;
         -m*n m*n m^2-n^2]
end

"""
    QI_layup(N::Integer)

Return vector with ply orientations [rad] for a symmetric quasi-isotropic
layup.

### Arguments
- N : Integer with the number of plies, must be multiple of 8.
"""
function QI_layup(N::Integer)
    if N % 8 == 0
        half = repeat([π/4, 0.0, -π/4, π/2], N÷8)
        return [half; reverse(half)]
    else
        throw(DomainError(N, "ply count must be multiple of eight."))
    end
end

"""
    CP_layup(N::Integer)

Return vector with ply orientations [rad] for a symmetric cross-ply layup.

### Arguments
- N : Integer with the number of plies, must be multiple of 8.
"""
function CP_layup(N::Integer)
    if N % 4 == 0
        half = repeat([0.0, π/2], N÷4)
        return [half; reverse(half)]
    else
        throw(DomainError(N, "ply count must be multiple of four."))
    end
end

"""
Type that represents a laminate.

### Fields

- plies::Vector{Ply}
    Layup.
"""
struct Laminate
    plies::Array{Ply}
end

function Laminate(mat::Material, layup::Vector{Float64}, t::Float64)
    plies = [Ply(mat, theta, t) for theta in layup]
    return Laminate(plies)
end


Base.iterate(L::Laminate, state=1) =
    state > length(L.plies) ? nothing : (L.plies[state], state+1)

Base.length(L::Laminate) = length(L.plies)

function Base.show(io::IO, L::Laminate)
    for ply in L
        println(ply)
    end
end



"""
    ply_interfaces(L::Laminate)

Return z-coordinates of ply interfaces of laminate `L`.
"""
function ply_interfaces(L::Laminate)
    z = zeros(length(L.plies) + 1)
    h = 0.0
    for (i, ply) in enumerate(L)
        h += ply.t
        z[i+1] = z[i] + ply.t
    end
    return z .- h/2
end

"""
    _repeated_ply_interfaces(L::Laminate)

Returns z-coordinates of the ply faces.
"""
function _repeated_ply_interfaces(L::Laminate)
    z = ply_interfaces(L)
    return [z[1]; repeat(z[2:end-1], inner=2); z[end]]
end

"""
    ABD_matrix(L::Laminate)

Calculates ABD-matrix for laminate `L`.
"""
function ABD_matrix(L::Laminate)
    A, B, D = zeros(3, 3), zeros(3, 3), zeros(3, 3)
    z = ply_interfaces(L)
    for (i, ply) in enumerate(L)
        A += stiffness_matrix(ply) * (z[i+1] - z[i])
        B += stiffness_matrix(ply) * (z[i+1]^2 - z[i]^2)/2
        D += stiffness_matrix(ply) * (z[i+1]^3 - z[i]^3)/3
    end
    return [A B; B D]
end

"""
    abd_matrix(L::Laminate)

Calculates abd-matrix for laminate `L`.
"""
function abd_matrix(L::Laminate)
    ABD = ABD_matrix(L)
    return inv(ABD)
end

"""
    engineering_constants(L::Laminate)

Calculates engineering constants of Laminate `L`.

### Argument

- L::Laminate

### Returns

- engcon::Dict{String, Float64}
    Dictionary with the following fields:
    - `Ex` : Young's modulus in x-direction.
    - `Ey` : Young's modulus in x-direction.
    - `Efx` : Flexural modulus in x-direction.
    - `Efy` : Flexural modulus in x-direction.
    - `Gxy` : Shear modulus.
    - `νxy` : Poisson modulus for loading in x-direction.
    - `νyx` : Poisson modulus for loading in y-direction.
"""
function engineering_constants(L::Laminate)
    abd = abd_matrix(L)
    z = ply_interfaces(L)
    H = z[end] - z[1]
    return Dict("Ex"=>1/(H*abd[1, 1]),
                "Ey"=>1/(H*abd[2, 2]),
                "Efx"=>12/(H^3*abd[4, 4]),
                "Efy"=>12/(H^3*abd[5, 5]),
                "Gxy"=>1/(H*abd[3, 3]),
                "νxy"=>-abd[2, 1]/abd[1, 1],
                "νyx"=>-abd[1, 2]/abd[2, 2])
end

"""
    Type that represents a loading case.

### Fields

- F::Vector{Float64}
    Force vector with `NaN` for the components that are not defined.
- d::Vector{Float64}
    Deformation vector with `NaN` for the components that are not defined.
- ΔT::Float64
    Temperature difference.

The force `F` and deformation vectors `d` have six elements:

1. Normal force `F[1]` and strain `d[1]` in X-direction
2. Normal force `F[2]` and strain `d[2]` in Y-direction
3. Shear force `F[3]` and shear strain `d[3]`
4. Bending moment `F[4]` and curavature `d[4]` along in X-direction
5. Bending moment `F[5]` and curavature `d[5]` along in Y-direction
6. Twisting moment `F[6]` and twisting curavature `d[6]`

A valid conditionrequires that only a force (or a moment) OR a strain
(or a curvature) is imposed in the given direction. As an example, in
case a normal force is applied in X-direction `F[1]`, the
corresponding strain `d[1]` cannot be provided and must be calculated.
The value for `d[1]` must therefore be equal to `NaN` in this case.

"""
struct Load
    F::Vector{Float64}
    d::Vector{Float64}
    ΔT::Float64
end

"""
    Load(F::Vector{Float64}, d::Vector{Float64})

Initializes a Load object with ΔT = 0.0.

### Arguments

- F::Vector{Float64}
    Force vector with `NaN` for the components that are not defined.
- d::Vector{Float64}
    Deformation vector with `NaN` for the components that are not defined.

"""
function Load(F::Vector{Float64}, d::Vector{Float64})
    Load(F, d, 0.0)
end

"""
    Load(load::Dict{String, Float64})

Initializes a Load from a dictionary.

### Arguments
- load::Dict
    Dictionary with the following keys:
    - `Fx` OR `ex`    : normal load OR strain in x-direction
    - `Fy` OR `ey`    : normal load OR strain in y-direction
    - `Fxy` OR `exy`  : shear load OR shear strain
    - `Mx` OR `kx`    : bending moment OR curvature in x-direction
    - `My` OR `ky`    : bending moment OR curvature in y-direction
    - `Mxy` OR `kxy`  : twisting moment OR curvature
    - `dT` (optional) : temperature difference
"""
function Load(load::Dict{String, Float64})
    labels = [("Fx", "ex"), ("Fy", "ey"), ("Fxy", "exy"),
              ("Mx", "kx"), ("My", "ky"), ("Mxy", "kxy")]
    F, d = NaN * ones(6), NaN * ones(6)
    for (i, label) in enumerate(labels)
        if haskey(load, label[1]) && !haskey(load, label[2])
            F[i] = load[label[1]]
        elseif !haskey(load, label[1]) && haskey(load, label[2])
            d[i] = load[label[2]]
        elseif haskey(load, label[1]) && haskey(load, label[2])
            throw(ErrorException(string("overdefined problem: both ",
                                        label[0],  " and ", label[1],
                                        " provided.")))
        elseif !haskey(load, label[1]) && !haskey(load, label[2])
            throw(ErrorException(string("overdefined problem: no ",
                                        label[0],  " or ", label[1],
                                        " provided.")))
        end
    end
    ΔT = haskey(load, "dT") ? load["dT"] : 0.0
    return Load(F, d, ΔT)
end

"""
    thermal_load(L::Laminate, ΔT::Float64)

Calculates thermal load and accompanying deformation vector.

### Arguments

- L::Laminate
    Laminate which is subjected to a temperature change ΔT.
- ΔT::Float64
    Temperature difference.

### Returns
- F_th::Vector{Float64}
    Fictive thermal load vector.
- d_th::Vector{Float64}
    Deformation due to thermal load.
"""
function thermal_load(L::Laminate, ΔT::Float64)
    z = ply_interfaces(L)
    N, M = zeros(3), zeros(3)
    for (i, ply) in enumerate(L)
        C, α = stiffness_matrix(ply), cte_vector(ply)
        N += C*α * ΔT * (z[i+1] - z[i])
        M += C*α * ΔT * (z[i+1]^2 - z[i]^2)/2
    end
    abd = abd_matrix(L)
    return [N; M], abd*[N; M]
end

"""
    thermal_load(L::Laminate, load::Load)

Calculates thermal load and accompanying deformation vector.

### Arguments

- L::Laminate
    Laminate which is subjected to a temperature change ΔT.
- load::Load
    Load object with ΔT.

### Returns
- F_th::Vector{Float64}
    Fictive thermal load vector.
- d_th::Vector{Float64}
    Deformation due to thermal load.

"""
function thermal_load(L::Laminate, load::Load)
    ΔT = load.ΔT
    thermal_load(L, ΔT)
end

"""
    apply_load(L::Laminate, load::Load)

"""
function apply_load(L::Laminate, load::Load)
    _, d_thermal = thermal_load(L, load)
    F, d = load.F, load.d
    iF, id = _load_masks(load)
    if all(iF)
        abd = abd_matrix(L)
        d = abd*F + d_thermal
    elseif all(id)
        ABD = ABD_matrix(L)
        F = ABD*(d - d_thermal)
    else
        abd = abd_matrix(L)
        abd_r = abd[:, iF][id, :]
        d_r = d[id] - abd_r*F[iF] - d_thermal[id]
        F[id] = abd[:, id][id, :] \ d_r
        d = abd*F + d_thermal
    end
    return F, d
end

"""
    _load_masks(load::Load)

Returns two masks with `1` for every non-NaN in `load.F` and `load.d`.

### Arguments

- load::Load
    Boundary conditions

### Returns
- F_mask::Vector{Bool}
    Mask for force vector.
- d_mask::Vector{Bool}
    Mask for deformation vector.

"""
function _load_masks(load::Load)
    return map(!isnan, load.F), map(!isnan, load.d)
end

"""
    is_valid(load::Load)

Returns true if boundary conditions in `load` are valid.

"""
function is_valid(load::Load)
    F, d = load_masks(load)
    return all(broadcast(⊻, F, d))  # ⊻ = XOR
end

"""
    torsion_shaft(T::Float64, R::Float64)

Returns Load object for a torsion shaft with radius `R` subjected to
torque `T`.

### Arguments

- T::Float64
    Torque.
- R::Float64
    Shaft radius.
"""
function torsion_shaft(T::Float64, R::Float64)
    F = [0, 0, T/(2π*R^2), 0, 0, 0]
    d = [NaN, NaN, NaN, NaN, NaN, NaN]
    return Load(F, d)
end

"""
    ply_strains(L::Laminate, load::Load)

Calculates strain on top and bottom surface of each ply in the ply coordinate
system.

### Arguments

- L::Laminate
    Laminate to which load is applied.
- load::Load
    Boundary conditions.

### Returns

- ε::Matrix{Float64}
    Strains in ply coordinate system. The returned matrix is of shape (2N, 3),
    with `N` the number of plies. For each ply the stress state at its top and
    bottom are returned. For the i-th ply:
    - column `i*2-1` holds the stress state at its top
    - column `i*2` holds the stress state at its bottom
    The three rows correspond to the three in-plane strains, i.e. the two
    normal strains and the shear strain.
- z::Vector{Float64}
    Vector with z-locations of the corresponding strains.
"""
function ply_strains(L::Laminate, load::Load)
    _, d = apply_load(L, load)
    z = _repeated_ply_interfaces(L)
    ε = broadcast(+, d[1:3], broadcast(*, d[4:6], z'))
    return ε, z
end

"""
    ply_stress(ply::Ply, ε::Vector{Float64}, ΔT::Float64=0.0)

Calculate ply stress as a function of strain and temperature change.

### Arguments

- ply::Ply
    Ply subjected to strain `ε`.
- ε::Vector{Float64}
    Strain vector.
- ΔT::Float64 (defaults to 0.0)
    Temperature difference.

### Returns

- σ::Vector{Float64}
    Stress vector.
"""
function ply_stress(ply::Ply, ε::Vector{Float64}, ΔT::Float64=0.0)
    C = stiffness_matrix(ply)
    α = cte_vector(ply)
    return C*(ε - α*ΔT)
end

"""
    ply_stresses(L::Laminate, load::Load)

Calculate ply stresses as function of boundary conditions.

### Arguments

- L::Laminate
    Laminate.
- load::Load
    Boundary conditions.

### Returns

- σ::Matrix{Float64}
    Stresses in ply coordinate system. The returned matrix is of shape (2N, 3),
    with `N` the number of plies. For each ply the stress state at its top and
    bottom are returned. For the i-th ply:
    - column `i*2-1` holds the stress state at its top
    - column `i*2` holds the stress state at its bottom
    The three rows correspond to the three in-plane stresses, i.e. the two
    normal stresses and the shear stress.
- z::Vector{Float64}
    Vector with z-locations of the corresponding strains.
"""
function ply_stresses(L::Laminate, load::Load)
    ε, z = ply_strains(L, load)
    ΔT = load.ΔT
    σ = zeros(size(ε))
    for (i, ply) in enumerate(L)
        σ[:, i*2-1] = ply_stress(ply, ε[:, i*2-1], ΔT)
        σ[:, i*2] = ply_stress(ply, ε[:, i*2], ΔT)
    end
    return σ, z
end

"""
    rotate_stress_to_matCS(σ::Vector{Float64}, θ::Float64)

Rotates a stress vector from ply to material coordinate system.

### Arguments

- σ::Vector{Float64}
    Stress vector in ply coordinate system.
- θ::Float64
    Rotation angle.

### Returns

- σₘ::Vector{Float64}
    Stress vector in material coordinate system.
"""
function rotate_stress_to_matCS(σ::Vector{Float64}, θ::Float64)
    T = transformation_matrix(θ)
    return T*σ
end

"""
    rotate_stress_to_matCS(σ::Vector{Float64}, ply::Ply)

Rotates a stress vector from ply to material coordinate system.

### Arguments

- σ::Vector{Float64}
    Stress vector in ply coordinate system.
- ply::Ply
    Ply object.

### Returns

- σₘ::Vector{Float64}
    Stress vector in material coordinate system.
"""
function rotate_stress_to_matCS(σ::Vector{Float64}, ply::Ply)
    return rotate_stress_to_matCS(σ, ply.θ)
end

"""
    rotate_stress_to_matCS(σ::Matrix{Float64}, L::Laminate)

Rotates a stress vector from ply to material coordinate system.

### Arguments

- σ::Matrix{Float64}
    Matrix with stress vectors on ply surfaces in ply CS.
- L::Laminate
    Laminate with layup information

### Returns

- σₘ::Matrix{Float64}
    Matrix with stress vectors in material coordinate system.
"""
function rotate_stress_to_matCS(σ::Matrix{Float64}, L::Laminate)
    if size(σ)[2] != 2*length(L)
        throw(ErrorException("Number of plies in `L` does not correspond to length of `σ`"))
    end
    σₘ = zeros(size(σ))
    for (i, ply) in enumerate(L)
        σₘ[:, i*2-1] = rotate_stress_to_matCS(σ[:, i*2-1], ply)
        σₘ[:, i*2] = rotate_stress_to_matCS(σ[:, i*2], ply)
    end
    return σₘ
end

"""
    pressure_vessel(ΔP::Float64, R::Float64)

Returns Load object for a pressure vessel with radius `R` and internal
pressure `ΔP`.

### Arguments

- ΔP::Float64
    Inside pressure minus external pressure.
- R::Float64
    Tank radius.
"""
function pressure_vessel(ΔP::Float64, R::Float64)
    F = [P*R/2, P*R, 0, 0, 0, 0]
    d = [NaN, NaN, NaN, NaN, NaN, NaN]
    return Load(F, d)
end
