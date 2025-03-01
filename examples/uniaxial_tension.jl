# Example: 24-ply QI laminate subjected to uniaxial load
#
# A quasi-isotropic fiber reinforced composite test specimen with a width of 25
# mm that is subjected to a tensile test. We will:
#
# 1. Determine the Young's modulus in X-direction
# 2. Check for failure in case a load of 50 kN is applied
using Plots
include("../clt.jl")

# First, we will need an instance of the Material class, which holds all
# relevant material properties. There are two ways to create a Material object.
# We can:
#
# 1. Import data from a json file
# 2. Provide a dictionary with the relevant data
#
# Here, we will import data from a json file.
fn = "materials/TC1200_UD.json"
TC1200 = Material(fn)

# Next we need to create a Laminate object. Again, there are two ways
# to do so. We can:
#
# 1. Provide a list of ply objects, which is demonstrated in
#    the next example
# 2. Provide a Material object, a list with fiber orientation angles,
#    and a ply thickness (equal for all plies).
#
# Here we use the second option, making use of the function `QI_layup()`
# to generate the orientation angles for a quasi-isotropic layup.
t = 0.15E-3
layup = QI_layup(24)
laminate = Laminate(TC1200, layup, t)

# The abd-matrix of the laminate is calculated as:
abd = abd_matrix(laminate)

# The engineering constants can be determined using:
eng_const = engineering_constants(laminate)

# This function yields a dictionary with the Young's and flexural
# moduli as well as the Poisson ratios of the laminate

# Next, a Load object is needed. There are three ways to create such an
# object. We can:
#
# 1. Import data from a json file
# 2. Provide a dictionary with the relevant data
# 3. Make use of a function for special loading cases
#
# Here, we will create a Load object from a dictionary. The dictionary
# should have the following keys:
# - Fx OR ex      : normal load or strain in x-direction
# - Fy OR ey      : normal load or strain in y-direction
# - Fxy OR exy    : shear load or shear strain
# - Mx OR kx      : bending moment or curvature in x-direction
# - My OR ky      : bending moment or curvature in y-direction
# - Mxy OR kxy    : twisting mometn or curvature
# - dT (optional) : temperature difference
F, width = 50E3, 25E-3
bc = Dict("Fx"=>F/width, "Fy"=>0.0, "Fxy"=>0.0,
          "Mx"=>0.0, "My"=>0.0, "Mxy"=>0.0, "dT"=>0.0)
uniaxial_tension = Load(bc)

# Now, we can determine the deformation of the laminate due to the
# applied load.
F, d = apply_load(laminate, uniaxial_tension)

# We can also calculate the stress distribution in the laminate.
σ_r, z = ply_stresses(laminate, uniaxial_tension)

# The resulting matrix σ_r holds the stresses distribution in ply
# coordinate system. It has `2N` columns and 3 rows. For the i-th ply,
# column `2i-1` holds the stress state at its top, while column `2i`
# holds the stress state at its bottom. The three rows correspond to
# the normal stresses in 1*, 2*, and the in-plane shear stress,
# respectively.

# The stresses can be rotated to the material coordinate system too.
σ = rotate_stress_to_matCS(σ_r, laminate)

# Plotting the stress in 1* direction in ply-coordinate system yields:
plot(σ[1, :]/1E6, z*1E3,
     xlabel="stress [MPa]",
     ylabel="z-coordinate [mm]",
     title="stress in 1*-direction",
     legend=false)
