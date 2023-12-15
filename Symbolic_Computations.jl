using Oscar, Primes

# This function computes the ED degree using msolve functionallity in Oscar.jl for a given tree and model
function EDdeg_symbolic(tree, model, projective=false, generic= false)

    A = AffineParametrization(tree, model)

    if projective
      A =  ProjectiveParametrization(A)
    end

    n, m = size(A)
    
    ps = Primes.primes(rand(10^2:10^4), rand(10^5:10^6))
    p = ps[rand(1:length(ps))]

    R, (x) = polynomial_ring(GF(p),["x$i" for i in 1:n+1])
    q = [prod(x[i]^A[i-1, j] for i in 2:n+1) for j in 1:m]

    # Distance function
    if generic
         λ, u = rand(1:10^7, m), rand(1:10^7, m)
        else
         λ, u = Vector{Int64}(ones(m)), rand(1:10^7, m) #λ = ones(m) should be changed for the proper values
    end
    g = sum(λ[i] * (q[i] - u[i])^2 for i in 1:length(q))

    # Construct system of critical equations
    F = [derivative(g, x[i]) for i in 2:n+1]
   
    f = convert(fpMPolyRingElem, prod(x) - 1)
    I = ideal(vcat(f,F))
    deg = Oscar.degree(ideal(groebner_basis_f4(I,eliminate=1))) # compute the degree of the ideal of critical equations after saturation
    return deg/cardinality(tree, model) # divide with the cardinality
end

# This function computes the cardinality of the fibers for a given tree and model
function cardinality(tree, model)

    intNodes = length(treeStructure(tree)[2]) + 1

    g = 1 

    if model in ["CFNmodel","K2model"]
        g = 2
    elseif model == "K3model"
        g = 4
    end

    return g^intNodes
end


