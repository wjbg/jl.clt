# Example: Flexural rigidity of a sandwich beam
#
# Calculate the flexural rigidity of sandwich beams. All beams have a
# width of 2 cm. Three different face sheet thicknesses are used,
# namely 0.25, 0.5 and 1.0 mm, while the core thickness is veried
# between 1 and 20 mm.
#
using Plots
include("../clt.jl")

# Core thicknesses, ply thickness and specimen width
t_core = 5E-3
t_skin = 0.25E-3

# First, we will need an instance of the Material class, which holds all
# relevant material properties. There are two ways to create a Material object.
# We can:
#
# 1. Import data from a json file
# 2. Provide a dictionary with the relevant data
#
# Here, we will provide a dictionary with the relevant data
CE = Material(Dict("name"=>"C/Epoxy",
                   "E1"=>120E9, "E2"=>9E9, "G12"=>4.8E9, "v12"=>0.29,
                   "alpha1"=>2E-6, "alpha2"=>40E-6,
                   "S1t"=>1800E6, "S1c"=>1200E6, "S2t"=>90E6, "S2c"=>100E6, "S6"=>105E6))

PU = Material(Dict("name"=>"PU foam",
                   "E1"=>250E6, "E2"=>250E6, "G12"=>95E6, "v12"=>0.32,
                   "alpha1"=>60E-6, "alpha2"=>60E-6,
                   "S1t"=>40E6, "S1c"=>50E6, "S2t"=>40E6, "S2c"=>50E6, "S6"=>25E6))

# Now we create two Ply objects, one for the plies in the C/epoxy and one for the core
cepoxy = Ply(CE, 0.0, t_skin)
core = Ply(PU, 0.0, t_core)

# Next we need to create a Laminate object. Again, there are two ways
# to do so. We can:
#
# 1. Provide a list of ply objects, which is demonstrated in
#    the next example
# 2. Provide a Material object, a list with fiber orientation angles,
#    and a ply thickness (equal for all plies).
#
# Here we use the first option to create three sandwich panels.
sandwich = [Laminate([cepoxy, core, cepoxy]),
            Laminate([cepoxy, cepoxy, core, cepoxy, cepoxy]),
            Laminate([cepoxy, cepoxy, cepoxy, core, cepoxy, cepoxy, cepoxy])]

# We can obtain the flexural modulus in x-direction as follows:
Efx = [engineering_constants(panel)["Efx"] for panel in sandwich]
