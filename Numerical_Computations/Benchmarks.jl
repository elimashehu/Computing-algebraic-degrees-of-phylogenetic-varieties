include("PhylogeneticData.jl");
include("Algebraic_Invariants_Functions.jl");

# Models: CFNmodel, JCmodel, K2Pmodel, K3Pmodel
# Trees: star3, star4, binary4, star5, binary5, tree5

tree = "star4"
model = "CFNmodel"


### EDdegrees ### Affine case
A = AffineParametrization(tree, model)

EDdeg_HC_A = EDdeg(A) #solving system of critical equations via solve function in HC.jl
tEDdeg_HC_A = @elapsed  EDdeg(A) # computing the time this computation take

EDdeg_M_A = monodromy_EDdeg(A) #solving system of critical equations via monodromy_solve function 
tEDdeg_M_A = @elapsed  monodromy_EDdeg(A)

# Generic EDdegrees # Affine case
gEDdeg_HC_A = EDdeg(A, true)
tgEDdeg_HC_A = @elapsed EDdeg(A, true)

gEDdeg_M_A = monodromy_EDdeg(A, true)
tgEDdeg_M_A = @elapsed monodromy_EDdeg(A, true)

### EDdegrees ### Projective case
P =  ProjectiveParametrization(A) #accepts two type of arguments: matrix or model and tree

EDdeg_HC_P = EDdeg(P)
tEDdeg_HC_P = @elapsed EDdeg(P)

EDdeg_M_P = monodromy_EDdeg(P)
tEDdeg_M_P = @elapsed monodromy_EDdeg(P)

# Generic EDdegrees # Projective case
gEDdeg_HC_P = EDdeg(P, true)
tEDdeg_HC_P = @elapsed EDdeg(P, true)

gEDdeg_M_P = monodromy_EDdeg(P,true)
tgEDdeg_M_P = @elapsed monodromy_EDdeg(P,true)

### MLdegrees ###
MLdeg_HC = MLdegree(tree, model) # using solve function in HC.jl
tMLdeg_HC = @elapsed MLdeg_HC

MLdeg_M = monodromy_MLdegree(tree, model) # using monodromy_solve function
tMLdeg_M = @elapsed MLdeg_M


# Perfroming this computations several times and saving the results in a text file

TREES = ["star3", "star4", "binary4", "star5", "binary5", "tree5"] # vector of strings defining the trees we consider
MODELS = ["CFNmodel", "JCmodel", "K2Pmodel", "K3Pmodel"]  # vector of strings defining the models we consider

file = string("results.txt") # we save the computed degrees in the file under the name results.txt

# we run the computations several times in this case 10 times
for i in 1:10
        for tree in TREES
                for model in MODELS
                        println(tree, " ", model, "\n")
                        ### EDdegrees ###
                        A = AffineParametrization(tree, model)

                        EDdeg_HC = EDdeg(A, model, false, true)
                        open(file, "a") do f
                                println(f, i, " ", tree, " ", model, " ", "ED", " ", EDdeg_HC)
                        end

                        MLdeg_M = monodromy_MLdegree(tree, model, true)
                        open(file, "a") do f
                                println(f, i, " ", tree, " ", model, " ", "ML", " ", MLdeg_M)
                        end
                end
        end
end
