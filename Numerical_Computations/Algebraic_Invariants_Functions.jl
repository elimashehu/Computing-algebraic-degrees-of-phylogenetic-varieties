using HomotopyContinuation, AbstractAlgebra, LinearAlgebra


function EDdeg(A, generic=false, feasibe=false, tolerance=1e-8)
    n, m = size(A)

    if generic
        λ, u = randn(Float64, m), rand(Float64, m)
    else
        λ, u = ones(m), rand(Float64, m) #λ = ones(m) should be changed for the proper values
    end

    # monomial parametrization
    @var θ[1:n]
    q=[prod(θ[i]^A[i,j] for i in 1:n) for j in 1:m]

    # Distance function
    g = sum(λ[i] * (q[i] - u[i])^2 for i in 1:length(q))

    # Construct system of critical equations
    F = System(differentiate(g,θ))

    # Find solutions
    R = solutions(HomotopyContinuation.solve(F, show_progress=true))
    nonSing = filter(r->all(abs.(r).>tolerance), R)

    # Extract the ceritified solutions along with biological meaningful ones 
    O = certify(F, nonSing)
    o = certificates(O)
    
    lll = map(o) do oᵢ
        if !is_certified(oᵢ)
            return false
        end
        𝕀 = certified_solution_interval(oᵢ)
        for i in 1:n
            if Base.in(0, HomotopyContinuation.IComplexF64(𝕀[i]))
                return false
            end
        end
        return true
    end

    nonSing = unique_points(nonSing)

    # Take unique solutions on the variety
    Q = System(q)
    qs = length(unique_points(map(s->Q(s), nonSing)))

    if feasibe
        solsProb, solsBio = feasibleSolutionsED(nonSing, model)
        qsProb = length(unique_points(map(s->Q(s), solsProb)))
        qsBio = length(unique_points(map(s->Q(s), solsBio)))
    else
        qsProb = qsBio = 0
    end

    return (qs, qsProb, qsBio, count(lll))
end

function monodromy_EDdeg(A, generic=false, tolerance=1e-8)
    n, m = size(A)

    if generic
        λ = randn(Float64, m)
    else
        λ = ones(m) #λ = ones(m) should be changed for the proper values
    end

    θ₀ = randn(n)
    u₀ = [prod(θ₀[i]^A[i,j] for i in 1:n) for j in 1:m]

    @var u[1:m]
    @var θ[1:n]

    q = [prod(θ[i]^A[i,j] for i in 1:n) for j in 1:m]

    # Distance function
    g = sum(λ[i] * (q[i] - u[i])^2 for i in 1:length(q))

    # Construct system of critical equations
    F = differentiate(g, θ)

    G = System(F, variables = θ, parameters = u)

    # Find solutions
    T = monodromy_solve(F, [θ₀], u₀, parameters = u)

    nonSing = filter(r -> all(abs.(r) .>tolerance), solutions(T))
    
    # Extract the ceritified solutions
    O = certify(G, nonSing, u₀)
    o = certificates(O)

    lll = map(o) do oᵢ
        if !is_certified(oᵢ)
            return false
        end
         𝕀 = certified_solution_interval(oᵢ)
            for i in 1:n
                if Base.in(0, HomotopyContinuation.IComplexF64(𝕀[i]))
                    return false
                end
            end
            return true
        end

    # Take unique solutions on the variety
    nonSing = unique_points(nonSing)
    Q = System(q)
    qs = map(s->Q(s),nonSing)

    return length(qs), count(lll)
end

# This function computes the ML degree using solve function for a given tree and model.
function MLdegree(tree, model, feasibe=false, tolerance=1e-8)
    # Compute the joint distribution at the leaves of the tree
    p = jointDistribution(treeStructure(tree), model)
    p = equivalentClasses(p)
    P = System(expand.(p))

    # List of parameters
    vars = HomotopyContinuation.variables(P)
    m = length(p); 
    n = length(vars);

    # Singular locus of the variety
    S2 = System(singularPoints(treeStructure(tree), model))

    # Take a random data point u
    u = rand(1:100, m)

    # List of parameters
    vars = HomotopyContinuation.variables(P)

    # System of critical equations.
    # New variables θ added to avoid rational expressions: θᵢpᵢ - 1 = 0
    @var θ[1:m]
    F1 = map(v -> sum(differentiate(p, v) .* u .* θ), vars)
    F2 = θ .* p .- ones(m)
    F = System(vcat(F1,F2))

    # Find all nonzero solutions
    R₀ = solutions(HomotopyContinuation.solve(F, show_progress=true))
    R = filter(r -> all(abs.(r) .>tolerance), R₀)

    # Take solutions with θᵢ, pᵢ ≠ 0
    G = System([prod(p)*prod(θ)])
    Rval = R[findall(x -> all(abs.(G(x)).>tolerance), R)]
    # Rval = R[findall(x -> all(abs.(convert(Array{ComplexF64}, expand.(G(x)))). >tolerance), R)]

    # Remove the θ
    Rₐ = map(x -> x[1:n], Rval)

    # Take the nonsingular points on the variety
    R₀ = Rₐ[findall(x -> all(abs.(S2(x)).>tolerance), Rₐ)]

    # Take the nonsingular real points
    R_real = filter(x -> all(isreal(round(y, digits=10)) for y in x), R₀)
    R_feasible = filter(r -> all(real.(r) .> -10^-10) && all(real.(r) .< 1 + 10^-10), R_real) 
    
    # Extract the ceritified solutions along with biological meaningful ones 
    H = differentiate(sum(u.*log.(p)), vars)
    O = certify(H, Rₐ)
    o = certificates(O)
    lll = map(o) do oᵢ
        if !is_certified(oᵢ)
            return false
        end
        𝕀 = certified_solution_interval(oᵢ)
        for i in 1:n
            if Base.in(0, HomotopyContinuation.IComplexF64(𝕀[i]))
                return false
            end
        end
        return true
    end

    # Take unique points on the variety
    points = unique_points(map(x -> P(x), (map(x -> x[1:n], R₀))))
    if feasibe  
        # points_feasible = unique_points(map(x -> P(x), (map(x -> x[1:n], R_feasible))))
        points_feasible = map(x -> P(x), (map(x -> x[1:n], R_feasible)))
        return length(points), count(lll), length(points_feasible)
    else
        return length(points), count(lll)
    end

end

# This function computes the ML degree using monodromy_solve function for a given tree and model.
function monodromy_MLdegree(tree, model, feasibe=false)
    # Compute the joint distribution at the leaves of the tree
    p = jointDistribution(treeStructure(tree), model)
    p = equivalentClasses(p)
    P = System(expand.(p))

    # Extract variables from the joint distribution system
    vars = HomotopyContinuation.variables(P)
    m = length(p); n = length(vars)
    r = min(n,m)

    @var u[1:m]

    # Compute the matrix A representing the differentiation of p with respect to the variables
    A = transpose(differentiate(p, vars))
    for i in 1:n
        A[i,:] = A[i,:]./p
    end

    # Generate random initial variables and find an intial solutions to the critical solutions associated to ML degree problem
    vars₀ = randn(Float64, length(vars))
    p₀ = HomotopyContinuation.evaluate(p, vars => vars₀)
    A₀ = convert(Matrix{ComplexF64}, expand.(Matrix{Expression}(HomotopyContinuation.evaluate(A, vars => vars₀))))
    u₁ = nullspace(A₀)[:, 1]
    
    # construct the system of critical equations
    G = System(A*u, variables = vars, parameters = u)
    T = monodromy_solve(G, [vars₀], u₁)
    sols = solutions(T)
    # nsols = length(sols)/fibers(tree, model)
    npoints, npointsProb, npointsBio = 0, 0, 0
    
    # Extract the ceritified solutions along with biological meaningful ones 
    if !isempty(sols)
        points = unique_points(map(x -> P(x), sols))
        npoints = length(points)

        G₀ = A * u₁
        O = certify(G₀, sols)
        o = certificates(O)
        zzz = map(o) do oᵢ
            if !is_certified(oᵢ)
                return false
            end
            𝕀 = certified_solution_interval(oᵢ)
            for i in 1:n
                if 0 in HomotopyContinuation.IComplexF64(𝕀[i])
                    return false
                end
            end
            return true
        end

        if feasibe
            solsProb, solsBio = feasibleSolutionsML(sols, model)
            # npointsProb = length(unique_points(map(x -> P(x), solsProb)))
            npointsProb = length(map(x -> P(x), solsProb))
            # npointsBio = length(unique_points(map(x -> P(x), solsBio)))
            npointsBio = length(map(x -> P(x), solsBio))
        end
    end

    if feasibe
        return npoints, npointsProb, npointsBio, count(zzz)
    else
        return npoints, count(zzz)
    end
end

# This function filters feasible solutions of Euclidean distance problem based on certain conditions for given solutions and a model.
function feasibleSolutionsED(sols, model)
    # Round and filter real solutions
    sols = [round.(sols[i]; digits=10) for i in 1:length(sols)]
    sols = sols[isreal.(sols)]
    
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