using HomotopyContinuation


#System of critical equations for ED problem
function ED_function(A, u, λ, n, m, θ)

    # monomial parametrization
    q = [prod(θ[i]^A[i,j] for i in 1:n) for j in 1:m]

    # Distance function
    g = sum(λ[i] * (q[i] - u[i])^2 for i in 1:length(q))

    # Construct system of critical equations
    return differentiate(g, θ), q  

end


#EDD function using solve function, output: number of certified solutions and the biological meaningful ones
function EDdeg(tree, model, projective=false, generic=false, feasibe=false, tolerance=1e-8)

    A, u, λ, n, m = generate_parameters(tree, model, projective, generic)

    @var θ[1:n] # Defining the variables
  
    # System of critical equations
    S, Q = ED_function(A, u, λ, n, m, θ)

    # Solve the system and extraxt the solutions
    R = solutions(HomotopyContinuation.solve(System(S), show_progress=true))

    nonSing = filter(r->all(abs.(r).>tolerance), R) # Nonsingular solutions
    nonSing = unique_points(nonSing) # Extract the solutions without repeation
 
    if !isempty(nonSing)
        # Certify the solutions
        O = certify(S, nonSing)
        # certificate_index(certificates(O)[1])#Check why this doesn't work
        cert_sols = ndistinct_certified(O)
        returnSolutions(model, System(Q), nonSing, cert_sols, feasibe)
    else
        return nothing
    end

end

#EDD function using monodromy, output: number of certified solutions and the biological meaningful ones
function monodromy_EDdeg(tree, model, projective=false, generic=false, feasibe=false, tolerance=1e-8, maxLoops=15)
    
    A, u, λ, n, m = generate_parameters(tree, model, projective, generic)

    # Define the variables
    @var θ[1:n] u[1:m]

    F, Q = ED_function(A, u, λ, n, m, θ) # Critical equations and parametrization

    # System of critical equations where u are parameters
    S = System(F, variables = θ, parameters = u)
    θ₁, u₁ = find_start_pair(S; max_tries = 100, atol = 0.0, rtol = 1e-12)

    # # Solve the system and extract the solutions
    T = solutions(monodromy_solve(F, [θ₁], u₁, parameters = u; max_loops_no_progress = maxLoops, group_actions = groupActions(model)))
   
    # Non Singular solutions with repeation
    nonSing = unique_points(filter(r -> all(abs.(r) .>tolerance), T))

    if !isempty(nonSing)  
        # certify the solutions
         O = certify(S, nonSing,  u₁)
        cert_sols = ndistinct_certified(O) 
        returnSolutions(model, System(Q), nonSing, cert_sols, feasibe)
    else
        return nothing
   end

end

function generate_parameters(tree, model, projective, generic)

    A = AffineParametrization(tree, model)

    if projective
        A =  ProjectiveParametrization(A)
    end

    n, m = size(A)

    if generic
        λ, u = randn(Float64, m), rand(Float64, m)
    else
        λ, u = ones(m), rand(Float64, m) #λ = ones(m) should be changed for the proper values
    end

    return A, u, λ, n, m 
end

function groupActions(model)
    if model == "CFNmodel"
        fCF(s) = -s
        return fCF
    elseif model == "K2Pmodel"
        fK2(s) = [(-1)^(i+1)*s[i] for i in 1:length(s)]
        return fK2
    elseif model == "K3Pmodel"
        exp1 = [-1, 1, -1]; f1(s) = [exp1[mod(i-1,3)+1]*s[i] for i in 1:length(s)]
        exp2 = [1, -1, -1]; f2(s) = [exp2[mod(i-1,3)+1]*s[i] for i in 1:length(s)]
        exp3 = [-1, -1, 1]; f3(s) = [exp3[mod(i-1,3)+1]*s[i] for i in 1:length(s)]
        return[f1,f2,f3]
    end
end

function f(tree, model)

    # Compute the joint distribution at the leaves of the tree
    p = jointDistribution(treeStructure(tree), model)
    p = equivalentClasses(p)
    P = System(HomotopyContinuation.expand.(p))
    
    # Extract variables from the joint distribution system
    vars = HomotopyContinuation.variables(P)
    m = length(p); 
    n = length(vars)

    return n, m, p, vars, P

end

# This function computes the ML degree using solve function for a given tree and model.
function MLdegree(tree, model, feasibe=false, tolerance=1e-8)
    
    n, m, p, vars, P = f(tree, model)

    # Take a random data point u
    u = rand(1:100, m)

    # System of critical equations.
    # New variables θ added to avoid rational expressions: θᵢpᵢ - 1 = 0
    @var θ[1:m]
    F1 = map(v -> sum(differentiate(p, v) .* u .* θ), vars)
    F2 = θ .* p .- ones(m)
    F =  System(vcat(F1,F2))

    # Find all nonzero solutions
    R₀ = solutions(HomotopyContinuation.solve(F, show_progress=true))
    R = filter(r -> all(abs.(r) .>tolerance), R₀)

    # Take solutions with θᵢ, pᵢ ≠ 0
    S = System([prod(p)*prod(θ)])
    Rval = R[findall(x -> all(abs.(S(x)).>tolerance), R)]

    # Remove the θ
    Rₐ = map(x -> x[1:n], Rval)

    # Singular locus of the variety
    S₁ = System(singularPoints(treeStructure(tree), model))
    # Take the nonsingular points on the variety
    R₁ = Rₐ[findall(x -> all(abs.(S₁(x)).>tolerance), Rₐ)]

    if !isempty(R₁)
       #certify the solutions
       H = differentiate(sum(u.*log.(p)), vars)
       O = certify(H, R₁)
       cert_sols = ndistinct_certified(O)

       returnSolutions(model, P, R₁, cert_sols, feasibe)
   else
     return nothing
   end 

end

# This function computes the ML degree using monodromy_solve function for a given tree and model.
function monodromy_MLdegree(tree, model, feasibe=false, maxLoops=5)

    n, m, p, vars, P = f(tree, model)
    
    # Compute the matrix B representing the differentiation of p with respect to the variables
    B = transpose(differentiate(p, vars))

       for i in 1:n
           B[i,:] = B[i,:]./p
       end
    
    # Generate random initial variables and find an intial solutions to the critical solutions associated to ML degree problem
    vars₀ = randn(Float64, length(vars))

    B₀ = convert(Matrix{ComplexF64}, HomotopyContinuation.expand.(Matrix{Expression}(HomotopyContinuation.evaluate(B, vars => vars₀))))
    u₁ = nullspace(B₀)[:, 1]
    
    @var u[1:m]
    # construct the system of critical equations
    S = System(B*u, variables = vars, parameters = u)
    T = monodromy_solve(S, [vars₀], u₁; max_loops_no_progress = maxLoops, group_actions = groupActions(model))
    sols = solutions(T)
    
    # Extract the ceritified solutions along with biological meaningful ones 
    if !isempty(sols)

        #certify the solutions 
        S₀ = B * u₁
        O = certify(S₀, sols)
        cert_sols = ndistinct_certified(O)
        
        returnSolutions(model, P, sols, cert_sols, feasibe)

    else

        return nothing

    end

end

# This function filters feasible solutions of Euclidean distance problem based on certain conditions for given solutions and a model.
function feasibleSolutionsED(sols, model)

    # Round and filter real solutions
    #sols=nonSing
    sols = [round.(sols[i]; digits=10) for i in 1:length(sols)]
    sols = sols[isreal.(sols)]
    print(length(sols))
    
    if length(sols) == 0
        return [], []
    else
         # Extract parameters and edges
        x = phyloParameters(model, 1)[1][1][1,:]
        param = length(unique(x[2:end]))
        edges = Int32(length(sols[1]) / param)
    
        solsProbs = falses(length(sols))
        solsBio = falses(length(sols))
    
        for i in 1:length(sols)
            S = reshape(Vector{Float64}(sols[i]), param, edges)
            M = zeros(ifelse(model == "CFNmodel", 1, 3), edges)
            
            # Construct M matrix based on the model
            if model == "CFNmodel"
                M[1,:] = (ones(edges) .- S[1,:]) ./ 2
            elseif model == "JCmodel"
                M[1,:] = M[2,:] = M[3,:] = (ones(edges) .- S[1,:]) ./ 4
            elseif model == "K2Pmodel"
                M[1,:] = M[2,:] = (ones(edges) .- S[1,:]) ./ 4   # b
                M[3,:] = (ones(edges) .- S[2,:]) ./ 2 .- M[1,:]  # a
            elseif model == "K3Pmodel"
                M[1,:] = (ones(edges) .+ S[1,:] .- S[2,:] .- S[3,:]) ./ 4   # a
                M[2,:] = (ones(edges) .- S[1,:] .+ S[2,:] .- S[3,:]) ./ 4   # b
                M[3,:] = (ones(edges) .- S[1,:] .- S[2,:] .+ S[3,:]) ./ 4   # c
            end
            
            M = [M; ones(1,edges) .- sum(M, dims=1)]
            
            # Check feasibility conditions
            if all(M .>= 0) && all(M .<= 1)
                solsProbs[i] = true
            end
            
            if all(M .>= 0) && diagonalLargeColumns(M)
                solsBio[i] = true
            end
        end
        return sols[solsProbs], sols[solsBio]
    end
end

# This function filters feasible solutions of Maximum likelihood problem based on certain conditions for given solutions and a model.
function feasibleSolutionsML(sols, model)
    
    # round the solutions
    sols = [round.(sols[i]; digits=10) for i in 1:length(sols)]
    
    # filter out complex solutions
    sols = sols[isreal.(sols)]

    x = phyloParameters(model, 1)[1][1][1,:] # extracting parameters from the model
    param = length(unique(x[2:length(x)])) # number of unique parameters

    edges = Int32(length(sols[1])/param) # number of edges in a tree
   
    # Initialize arrays to store feasibility results for probability and biological constraints.
    solsProbs = falses(length(sols))
    solsBio = falses(length(sols))

     # Iterate over each solution to check feasibility.
    for i in 1:length(sols)
        S = reshape(Vector{Float64}(sols[i]), param, edges)
       
        #Construct matrix M depending on the model
        if model == "CFNmodel"
            M = [S; reshape(ones(edges) - S[1,:], 1, edges)]
        elseif model == "JCmodel"
            M = [S; reshape(ones(edges) - 3*S[1,:], 1, edges)]
        elseif model == "K2Pmodel"
            M = [S; reshape(ones(edges) - S[1,:] - 2*S[2,:], 1, edges)]
        elseif model == "K3Pmodel"
            M = [S; reshape(ones(edges) - S[1,:] - S[2,:] - S[3,:], 1, edges)]
        end

         # Check if all elements of M satisfy probability constraints.
        if all(M .>= 0) & all(M .<= 1)
            solsProbs[i] = true
        end

         # Check if all elements of M satisfy diagonal and large column constraints.
        if all(M .>= 0) & diagonalLargeColumns(M)
            solsBio[i] = true
        end
    end

    # Return real solutions different from zero that satisfy probability and biological constraints.
    return(sols[solsProbs], sols[solsBio])

end

# Define a helper function to check if diagonals are larger than columns
function diagonalLargeColumns(M)

 rows = size(M)[1]
   
    return(all(transpose(M[rows,:]) .>= M[1:rows-1,:]))

end

# Evaluate a list of points into a system and take the unique ones
function unique_evaluate(Sys, points)
    if length(points) != 0
        uniquePoints = map(x -> Sys(x), unique_points(points))
        uniquePoints = unique_points([convert(Vector{ComplexF64}, uniquePoints[i]) for i in 1:length(uniquePoints)])
        return uniquePoints
    else
        return []
    end
end

function returnSolutions(model, Sys, nonSing, cert_sols, feasibe)
    # Take unique solutions on the variety
    npoints = length(unique_evaluate(Sys, nonSing))
    #npoints = length(unique_points(map(s->Sys(s), nonSing))) # Use certified sols and not nonsing one?
    
    # Extract the ceritified and the biologically meaningful solutions
    if feasibe
        solsProb, solsBio = feasibleSolutionsED(nonSing, model)
        npointsProb = length(unique_evaluate(Sys, solsProb))
        npointsBio = length(unique_evaluate(Sys, solsBio))
        returnData = (cert_sols, npoints, npointsProb, npointsBio)
    else
       returnData = (cert_sols, npoints)
    end
end
