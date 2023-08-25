using HomotopyContinuation, AbstractAlgebra, LinearAlgebra, IterTools
using Combinatorics

# Models: CFNmodel, JCmodel, K2Pmodel, K3Pmodel
# Trees: star3, star4, binary4, star5, binary5, tree5

# Define matrices of the monomial parametrization of the affine Varieties for different models and trees
M_T = Dict("CFNmodel" =>
Dict(
"star3" => [0 1 1; 1 0 1; 1 1 0],
"star4" => [0 0 0 1 1 1 1; 0 1 1 0 0 1 1; 1 0 1 0 1 0 1; 1 1 0 1 0 0 1],
"binary4" => [0 1 1 1 1 0 0; 0 0 0 1 1 1 1; 0 1 1 0 0 1 1; 1 0 1 0 1 0 1; 1 1 0 1 0 0 1],
"star5" => [0 0 0 0 0 0 0 1 1 1 1 1 1 1 1; 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1; 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1; 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1; 1 1 0 1 0 0 1 1 0 0 1 0 1 1 0],
"binary5" => [0 0 0 1 1 1 1 1 1 1 1 0 0 0 0; 0 1 1 1 1 0 0 1 1 0 0 0 0 1 1; 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1; 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1; 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1; 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1; 1 1 0 1 0 0 1 1 0 0 1 0 1 1 0],
"tree5" => [0 0 0 1 1 1 1 1 1 1 1 0 0 0 0; 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1; 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1; 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1; 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1; 1 1 0 1 0 0 1 1 0 0 1 0 1 1 0]),
"JCmodel" =>
Dict(
"star3" => [0 1 1 1; 1 0 1 1; 1 1 0 1],
"star4" => [0 0 0 0 1 1 1 1 1 1 1; 0 1 1 1 0 0 0 1 1 1 1; 1 0 1 1 0 1 1 0 1 0 1; 1 1 0 1 1 0 1 0 1 1 0],
"binary4" => [0 1 1 1 1 1 1 0 0 1 1 1; 0 0 0 0 1 1 1 1 1 1 1 1; 0 1 1 1 0 0 0 1 1 1 1 1; 1 0 1 1 0 1 1 0 1 0 1 1; 1 1 0 1 1 0 1 0 1 1 1 0],
"star5" => [0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1; 0 1 1 1 0 0 0 1 1 1 1 0 0 0 1 1 1 1 0 0 1 1 1 0 0 1; 1 0 1 1 0 1 1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 1 0 1 0; 1 1 0 1 1 0 1 0 1 1 0 1 0 1 0 1 1 0 0 1 1 0 1 1 0 0],
"binary5" => [0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 0 0 1 1 1 1 1 1 0 0 1 1 1 0 0 1 1 1 1 1 1 1 1 1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1; 0 1 1 1 0 0 0 1 1 1 1 1 0 0 0 1 1 1 1 1 0 0 1 1 1 0 0 0 1 1 1 1 1; 1 0 1 1 0 1 1 0 1 0 1 1 0 1 1 0 1 0 1 1 0 1 0 1 1 0 1 1 0 1 1 0 1; 1 1 0 1 1 0 1 0 1 1 1 0 1 0 1 0 1 1 1 0 0 1 1 0 1 1 1 0 1 1 0 0 1],
"tree5" => [0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1; 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1; 0 1 1 1 0 0 0 1 1 1 1 0 0 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1 1 1; 1 0 1 1 0 1 1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 1 0 1 1 0 1 1 0; 1 1 0 1 1 0 1 0 1 1 0 1 0 1 0 1 1 0 0 1 1 0 1 1 1 0 1 1 0 0]),
"K2Pmodel" =>
Dict(
"star3" => [0 0 1 1 1 0 0 0 0; 0 0 0 0 0 1 1 1 1; 1 0 0 1 0 0 1 0 0; 0 1 0 0 1 0 0 1 1; 1 0 1 0 0 0 0 0 1; 0 1 0 0 1 1 1 0 0],
"star4" => [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1; 0 0 1 1 1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0; 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1; 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 1; 0 1 0 0 1 0 0 1 1 0 0 1 0 0 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 0 0; 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 1 0 0 0 0 1 0 0 1 0 0 1 0 1 0; 0 1 0 0 1 1 1 0 0 0 0 1 0 0 1 1 1 0 0 1 1 0 0 1 1 0 0 0 0 1 0 0],
"binary4" => [0 0 1 1 1 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1; 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 0 0 1 1 1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0; 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1; 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 1 0; 0 1 0 0 1 0 0 1 1 0 0 1 0 0 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 0 0 1; 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 1 0 0 0 0 1 0 0 1 0 0 1 0 1 0 0; 0 1 0 0 1 1 1 0 0 0 0 1 0 0 1 1 1 0 0 1 1 0 0 1 1 0 0 0 0 1 0 0 1],
"star5" => [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 0 0 1 1 1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 1 1; 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0; 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 1 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0 1 0 0 1 0 1 0 0 1 0 0 1 0 0 0 1 0 1; 0 1 0 0 1 0 0 1 1 0 0 1 0 0 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 0 0 0 0 1 0 0 1 0 0 1 1 0 0 1 0 0 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 0 0 0 0 1 1 0 0 1 1 0 0 1 0 0 0 0 1 1 0 0 1 1 0 0 1 0 0 0 0 1 0 0 1 0 0 1 1 0 0 0 0; 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 1 0 0 0 0 1 0 0 1 0 0 1 0 1 0 1 0 0 0 1 0 0 0 1 0 0 1 0 1 0 0 0 0 0 1 0 0 1 0 0 0 0 1 1 0 0 0 1 0 0 0 1 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0 0 1 1 0 0 0 1 0 1 0 1 0 0 0 0 0 1 1 0 0 1; 0 1 0 0 1 1 1 0 0 0 0 1 0 0 1 1 1 0 0 1 1 0 0 1 1 0 0 0 0 1 0 0 0 0 1 0 0 1 1 1 0 0 0 0 1 0 0 1 1 1 0 0 1 1 0 0 1 1 0 0 0 0 1 0 0 1 1 0 0 1 1 0 0 0 0 1 0 0 1 1 0 0 1 1 0 0 0 0 1 0 0 0 0 1 0 0 1 1 1 0 0 0 0 0 0],
"binary5" => [],
"tree5" => []),
"K3Pmodel" =>
Dict(
"star3" => [0 0 0 1 1 1 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1; 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0; 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0; 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1; 1 0 0 1 0 0 0 0 0 0 1 0 0 1 0; 0 1 0 0 0 0 1 1 0 0 0 0 1 0 0; 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0],
"star4" => [],
"binary4" => [],
"star5" => [],
"binary5" => [],
"tree5" => [])
)

# Given a tree and a model as inputs, return the corresponding A matrix
function AffineParametrization(tree, model)
    keys(M_T)
     # Check if the model is present in the defined data
    if !(model in keys(M_T))
        return "Model not defined"
    end
    # Check if the tree is present in the defined data
    if !(tree in keys(M_T[model]))
        return "Tree not defined"
    end
    # Check if the data associated to the model and tree is empty
    if isempty(M_T[model][tree])
        return "No data found"
    end
    # If all checks pass, return the matrix
    return M_T[model][tree]
end
# The following two functions do accept two type of arguments
# But both of them augment the matrix A with additional row and column to create a projective map
# Vertically concatenate A with a column of zeros
# Horizontally concatenate A with a row of ones
function ProjectiveParametrization(A)
    (n, m) = size(A)
    Matrix{Int64}(vcat(hcat(A, zeros(n)), transpose(ones(m+1))))
end
function ProjectiveParametrization(tree, model)
    # gives a parametrising map of the projective variety
    A = AffineParametrization(tree, model)
    (n, m) = size(A)
    vcat(hcat(A, zeros(n)), ones(1, m+1))
end
function CFmatrices(n)
    @var a[1:n]
    M = repeat([Array{Expression}(undef, 2, 2)], n)
    M = [[1-a[i] a[i]; a[i] 1-a[i]] for i in 1:n]
end

function GBmatrices(model, n)
    @var a[1:n];
    if model == "K2Pmodel"
        @var b[1:n]
        c = b
    elseif model == "JCmodel"
        b = a; c = a;
    else
        @var b[1:n]; @var c[1:n];
    end

    M = repeat([Array{Expression}(undef, 4, 4)], n)
    for i in 1:n
        M[i] = [1-(a[i]+b[i]+c[i]) a[i] b[i] c[i]; a[i] 1-(a[i]+b[i]+c[i]) c[i] b[i]; b[i] c[i] 1-(a[i]+b[i]+c[i]) a[i]; c[i] b[i] a[i] 1-(a[i]+b[i]+c[i])]
    end
    M
end

function phyloParameters(model, nEdges)
    if model == "CFNmodel"
        CFmatrices(nEdges), 2
    else
        GBmatrices(model, nEdges), 4
    end
end

function jointDistribution(tree, model)
    leaves = tree[1]
    nL = length(leaves)
    intP = tree[2]

    nEdges = nL + length(intP)
    M, param = phyloParameters(model, nEdges)

    p = reshape(Array{Expression}(undef, param^nL), Tuple(ones(Int64, nL) * param))
    indices = CartesianIndices(p)

    for i in 1:length(indices)
        total_sum = zero(Expression)

        if nL == 3
            sum_expr = sum(M[j][:, indices[i][j]] for j in 1:nL)
            total_sum = 1 / param * prod(sum_expr)
        elseif nL == 4
            if nEdges == 4
                sum_expr = sum(M[j][:, indices[i][j]] for j in 1:nL)
                total_sum = 1 / param * prod(sum_expr)
            else
                mat_expr = prod(M[j][:, indices[i][j]] for j in 1:2)
                sum_expr = sum(M[j][:, indices[i][j]] for j in 3:nL)
                total_sum = 1 / param * transpose(mat_expr) * M[nL] * sum_expr
            end
        elseif nL == 5
            if nEdges == 5
                sum_expr = sum(M[j][:, indices[i][j]] for j in 1:nL)
                total_sum = 1 / param * prod(sum_expr)
            elseif nEdges == 6
                mat_expr = prod(M[j][:, indices[i][j]] for j in 1:2)
                sum_expr = sum(M[j][:, indices[i][j]] for j in 3:nL)
                total_sum = 1 / param * transpose(mat_expr) * M[nL] * sum_expr
            elseif nEdges == 7
                inner_sum = zero(Expression)
                for y1 in 1:param, y2 in 1:param, y3 in 1:param
                    inner_sum += M[1][y1, indices[i][1]] * M[2][y1, indices[i][2]] * M[6][y1, y2] *
                                 M[7][y2, y3] * M[3][y2, indices[i][3]] * M[4][y3, indices[i][4]] *
                                 M[5][y3, indices[i][5]]
                end
                total_sum = 1 / param * inner_sum
            end
        end

        p[i] = total_sum
    end

    p
end



function equivalentClasses(p)
    p = reshape(p, length(p))
    eq = unique(p)

    for i in 1:length(eq)
        eq[i] = eq[i] * length(findall(x -> x == eq[i], p))
    end
    eq
end


function singularPoints(tree, model)
    leaves = tree[1]
    nL = length(leaves)
    intP = tree[2]

    nEdges = nL + length(intP)
    M, param = phyloParameters(model, nEdges)

    sing = []  # Preallocate the sing vector

    if model == "CFNmodel"
        for i in 1:length(M)
            λ₁ = M[i][1, 1] - M[i][1, 2]
            push!(sing, λ₁)  # Use push! to add elements to the vector
        end
    else
        for i in 1:length(M)
            λ₁ = M[i][1, 1] + M[i][1, 2] - M[i][1, 3] - M[i][1, 4]
            λ₂ = M[i][1, 1] - M[i][1, 2] + M[i][1, 3] - M[i][1, 4]
            λ₃ = M[i][1, 1] - M[i][1, 2] - M[i][1, 3] + M[i][1, 4]
            append!(sing, [λ₁, λ₂, λ₃])  # Use append! to concatenate arrays
        end
    end

    return unique(sing)  # Return the unique singular points
end


# Define the tree structures using a dictionary
tree_structures = Dict(
    "star3" => [[1, 2, 3], []],
    "star4" => [[1, 2, 3, 4], []],
    "binary4" => [[1, 2, 3, 4], [[1, 2]]],
    "star5" => [[1, 2, 3, 4, 5], []],
    "binary5" => [[1, 2, 3, 4, 5], [[1, 2], [1, 2, 3]]],
    "tree5" => [[1, 2, 3, 4, 5], [[1, 2]]]
)

# Function to get the tree structure based on the tree name
function treeStructure(tree)
    if haskey(tree_structures, tree)
        return tree_structures[tree]
    else
        return "Tree not defined"
    end
end


