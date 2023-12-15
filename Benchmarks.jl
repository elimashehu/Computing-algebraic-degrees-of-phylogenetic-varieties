include("PhylogeneticData.jl");
include("Algebraic_Invariants_functions.jl");
include("Symbolic_Computations.jl");

# Models: CFNmodel, JCmodel, K2Pmodel, K3Pmodel
# Trees: star3, star4, binary4, star5, binary5, tree5
tree = "star3"
model = "CFNmodel"


### EDdegrees ### Affine case

# Symbolic Computations in Oscar.jl
Edd_msolve = EDdeg_symbolic(tree, model)
gEdd_msolve = EDdeg_symbolic(tree, model, false, true)

# Numerical Computations in HomotopyContinuation.jl
EDdeg_HC_A = EDdeg(tree,model) #solving system of critical equations via solve function in HC.jl

# Generic EDdegree
gEDdeg_HC_A = EDdeg(tree, model, false, true)


### EDdegrees ### Projective case

# Symbolic Computations in Oscar.jl
Edd_msolve = EDdeg_symbolic(tree, model, true)
gEdd_msolve = EDdeg_symbolic(tree, model, true, true)

# Numerical Computations in HomotopyContinuation.jl
EDdeg_HC_P = EDdeg(tree,model,true)

# Generic EDdegree
gEDdeg_M_P = monodromy_EDdeg(tree,model,true,true) #solving system of critical equations via monodromy_solve function 


### MLdegrees ###
MLdeg_HC = MLdegree(tree, model) # using solve function in HC.jl

MLdeg_M = monodromy_MLdegree(tree, model) # using monodromy_solve function
