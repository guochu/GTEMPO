abstract type ExponentialExpansionAlgorithm end
abstract type AbstractPronyExpansion <: ExponentialExpansionAlgorithm end

# hankel expansion
struct PronyExpansion <: AbstractPronyExpansion
    n::Int 
    stepsize::Int
    tol::Float64
    verbosity::Int
end
PronyExpansion(; n::Int=10, stepsize::Int=1, tol::Real = 1.0e-8, verbosity::Int=1) = PronyExpansion(n, stepsize, convert(Float64, tol), verbosity)

struct DeterminedPronyExpansion <: AbstractPronyExpansion
    n::Int 
    stepsize::Int
    tol::Float64
    verbosity::Int
end
DeterminedPronyExpansion(; n::Int=10, stepsize::Int=1, tol::Real = 1.0e-8, verbosity::Int=1) = DeterminedPronyExpansion(n, stepsize, convert(Float64, tol), verbosity)


function prony(x::Vector, p::Int)
    n = length(x)
    @assert p <= n÷2 "p can not exceed length(x)/2"

    # find the roots of characteristic polynomial
    A = zeros(typeof(x[1]), p, p)
    for i = 1:p, j = 1:p
        A[i,j] = x[p+i-j]
    end
    a = -A\x[p+1:2p]
    pushfirst!(a,1)
    z = roots(Polynomial(reverse(a)))

    # find the coefficient
    A = zeros(typeof(z[1]), p, p)
    for i = 1:p, j = 1:p
        A[i,j] = z[j]^i
    end
    α = A\x[1:p]

    # # compute the error
    # y = zeros(typeof(α[1]), n)
    # for i = 1:n, k = 1:p
    #     y[i] += α[k]*z[k]^i
    # end
    # E = norm(x-y)
        
    α, z
end

# least square version
function lsq_prony(x::Vector, p::Int)
    n = length(x)
    @assert p <= n÷2 "p can not exceed length(x)/2"
    
    # find the roots of characteristic polynomial via least square method
    A = zeros(typeof(x[1]), n-p , p)
    for i = 1:n-p, j = 1:p
        A[i,j] = x[p+i-j]
    end
    a = -A\x[p+1:n]
    pushfirst!(a,1)
    z = roots(Polynomial(reverse(a)))

    # find the coefficient via least square method
    A = zeros(typeof(z[1]), n, p)
    for i = 1:n, j = 1:p
        A[i,j] = z[j]^i
    end
    α = A\x

    # # compute the error
    # y = zeros(typeof(α[1]), n)
    # for i = 1:n, k = 1:p
    #     y[i] += α[k]*z[k]^i
    # end
    # E = norm(x-y)
        
    α, z
end

exponential_expansion_n(f::Vector, p::Int, alg::PronyExpansion) = lsq_prony(f, p)
exponential_expansion_n(f::Vector, p::Int, alg::DeterminedPronyExpansion) = prony(f, p)

# function exponential_expansion_n(f::Vector, p::Int, alg::AbstractPronyExpansion)
#     α, z, E = _exponential_expansion_n(f, p, alg)
#     for i in 1:length(α)
#         α[i] /= z[i]
#     end
#     return α, z, E
# end

function exponential_expansion(f::Vector{<:Number}, alg::AbstractPronyExpansion)
    (length(f) > 1) || throw(ArgumentError("length of data should be larger than 1"))
    xs, lambdas = _exponential_expansion_impl(f, alg)
    if alg.stepsize != 1
        expansion_changestepsize!(xs, lambdas, alg.stepsize)
    end
    if alg.verbosity > 2
        println("Prony coefs: ", xs)
        println("Prony roots: ", lambdas)
    end
    return xs, lambdas 
end

function expansion_changestepsize!(xs::Vector, lambdas::Vector, stepsize::Int)
    α = 1.0/stepsize
    for i in 1:length(lambdas)
        xs[i] *= (lambdas[i])^(1-α)
        lambdas[i] = (lambdas[i])^(α)
    end    
    return xs, lambdas
end

function _exponential_expansion_impl(f::Vector{<:Number}, alg::AbstractPronyExpansion)
    L = length(f)
    tol = alg.tol
    verbosity = alg.verbosity
    maxiter = alg.n
    nitr = min(maxiter, L)
    for n in 1:nitr
        xs, lambdas = exponential_expansion_n(f, n, alg)
        err = expansion_error(f, xs, lambdas)
        if err <= tol
            (verbosity > 1) && println("PronyExpansion converged in $n iterations, error is $err")
            # println(xs, " ", lambdas)
            return xs, lambdas
        end
        if n >= min(L-n+1, nitr)
            (verbosity > 0) && @warn "can not find a good approximation with L=$(L), n=$(alg.n), tol=$(tol), return with error $err"
            # println(xs, " ", lambdas)
            return xs, lambdas
        end
    end
    error("can not be here")
    # @warn "can not find a good approximation with n=$(maxiter), tol=$(tol), try increase L, or decrease tol"
end

function _predict(x, p)
    @assert length(p) % 2 == 0
    n = div(length(p), 2)
    L = length(x)
    T = eltype(p)
    r = zeros(T, L)
    for i in 1:L
        xi = x[i]
        @assert xi == i
        tmp = zero(T)
        for j in 1:n
            tmp += p[j] * p[n+j]^xi
        end
        r[i] = tmp
    end
    return r
end

function expansion_error(f::Vector{<:Number}, p::Vector{<:Number})
    T = eltype(f)
    xdata = [convert(T, i) for i in 1:length(f)]
    f_pred = _predict(xdata, p)
    return norm(f_pred - f)
end
expansion_error(f::Vector{<:Number}, coeffs::Vector{<:Number}, alphas::Vector{<:Number}) = expansion_error(f, vcat(coeffs, alphas))

# # least square expansion
# struct LsqExpansion <: ExponentialExpansionAlgorithm 
#     tol::Float64
#     verbosity::Int
# end
# LsqExpansion(; tol::Real = 1.0e-8, verbosity::Int=1) = LsqExpansion(convert(Float64, tol), verbosity)


# function lsq_expansion_n(f::Vector{<:Real}, n::Int, coeffs::Vector{<:Real}, alphas::Vector{<:Real})
#     @assert n == length(coeffs) == length(alphas)
#     T = eltype(f)
#     xdata = [convert(T, i) for i in 1:length(f)]
#     p0 = vcat(coeffs, alphas)
#     fit = curve_fit(_predict, xdata, f, p0, autodiff=:forwarddiff)
#     # println("converged? ", fit.converged)
#     p = fit.param
#     err = norm(_predict(xdata, p) - f)
#     return p[1:n], p[n+1:end], err
# end
# function exponential_expansion(f::Vector{<:Number}, alg::LsqExpansion)
#     L = length(f)
#     results = []
#     errs = Float64[]
#     tol = alg.tol
#     verbosity = alg.verbosity
#     for n in 1:L
#         if n == 1
#             xs, lambdas, err = lsq_expansion_n(f, n, [0.5], [0.5])
#         else
#             _xs, _lambdas = results[end]
#             xs, lambdas, err = lsq_expansion_n(f, n, vcat(_xs, 0.5), vcat(_lambdas, 0.5))
#         end
#         if err <= tol
#             (verbosity > 1) && println("converged in $n iterations, error is $err")
#             return xs, lambdas
#         else
#             push!(results, (xs, lambdas))
#             push!(errs, err)
#         end
#         if n == L
#             (verbosity > 0) && @warn "can not converge to $tol with size $L, try increase L, or decrease tol"
#             return results[argmin(errs)]
#         end
#     end
#     error("can not find a good approximation")
# end

"""
    exponential_expansion(f::Vector{<:Number}; alg::ExponentialExpansionAlgorithm=PronyExpansion())

Return a list of αᵢ and βᵢ which satisfy:
f(x) = ∑ᵢ αᵢ × (βᵢ)ˣ, for 1 ≤ x ≤ N
"""
exponential_expansion(f::Vector{<:Number}; alg::ExponentialExpansionAlgorithm=PronyExpansion()) = exponential_expansion(f, alg)
exponential_expansion(f, L::Int, alg::ExponentialExpansionAlgorithm) = exponential_expansion([f(k) for k in 1:L], alg)
exponential_expansion(f, L::Int; alg::ExponentialExpansionAlgorithm=PronyExpansion()) = exponential_expansion(f, L, alg)
