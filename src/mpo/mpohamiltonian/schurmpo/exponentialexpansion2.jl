
# hankel expansion
struct PronyExpansion2 <: AbstractPronyExpansion
    n::Int 
    atol::Float64
    rtol::Float64
    verbosity::Int
end
PronyExpansion2(; n::Int=10, atol::Real = 1.0e-8, rtol::Real = 1.0e-3, verbosity::Int=1) = 
            PronyExpansion2(n, convert(Float64, atol), convert(Float64, rtol), verbosity)
Base.similar(alg::PronyExpansion2; n::Int=alg.n, atol::Real=alg.atol, rtol::Real=alg.rtol, verbosity::Int=alg.verbosity) = 
            PronyExpansion2(n=n, atol=atol, rtol=rtol, verbosity=verbosity)


struct LsqExpansion2 <: ExponentialExpansionAlgorithm
    n::Int
    atol::Float64
    rtol::Float64
    verbosity::Int
end
LsqExpansion2(; n::Int=10, atol::Real = 1.0e-8, rtol::Real = 1.0e-3, verbosity::Int=1) = 
            LsqExpansion2(n, convert(Float64, atol), convert(Float64, rtol), verbosity)
Base.similar(alg::LsqExpansion2; n::Int=alg.n, atol::Real=alg.atol, rtol::Real=alg.rtol, verbosity::Int=alg.verbosity) = 
            LsqExpansion2(n=n, atol=atol, rtol=rtol, verbosity=verbosity)


const ExponentialExpansionAlgorithm2 = Union{PronyExpansion2, LsqExpansion2}



function exponential_expansion(f::Vector{<:Number}, alg::ExponentialExpansionAlgorithm2)
	r_atol = norm(f)*alg.rtol
    if r_atol < alg.atol
        (alg.verbosity >= 1) && println("Using atol of $r_atol according to rtol")
        alg = similar(alg, atol=r_atol)
    end

    
    idx = first_period(real(f))
    steps = unique([[round(Int, idx*i) for i in [0.2, 0.3, 0.35, 0.4, 0.45]]; 1])
    filter!(x->(x>0), steps)
    xs, lambdas = _exponential_expansion(f, alg, steps=steps)
    xs, lambdas = cut(f, xs, lambdas, alg)
    if alg.verbosity >= 2
        println("Prony coefs: ", xs)
        println("Prony roots: ", lambdas)
    end
    return xs, lambdas 
end



# prony
function _exponential_expansion(f::Vector{<:Number}, alg::PronyExpansion2; steps::Vector{Int})
    L = length(f)
    atol = alg.atol
    nitr = min(alg.n, L)
    for n in 1:nitr
        (xs, lambdas), err = exponential_expansion_n(f, n, alg, steps=steps)
        if err <= atol
            (alg.verbosity >= 1) && println("PronyExpansion2 converged in $n iterations, error is $err")
            return xs, lambdas
        end
        if n >= min(L-n+1, nitr)
            @warn "can not find a good approximation with L=$(L), n=$(alg.n), atol=$(atol), return with error $err"
            return xs, lambdas
        end
    end
    error("can not be here")
end

function exponential_expansion_n(f::Vector, p::Int, alg::PronyExpansion2; steps::Vector{Int})
    errs = []
    coeffs = []
    for step in steps
        f2 = f[1:step:end]
        (p > length(f2)÷2) && continue
        
        xs, lambdas = lsq_prony(f2, p)
        xs′ = @. xs * lambdas ^(1-1/step)
        lambdas′ = @. lambdas ^ (1/step)
        err = expansion_error(f, xs′, lambdas′)

        push!(errs, err)
        push!(coeffs, (xs′, lambdas′))
    end

    _, idx = findmin(identity, errs)
    return coeffs[idx], errs[idx]
end



# lsq
_exponential_expansion(f::Vector{<:Real}, alg::LsqExpansion2; steps::Vector{Int}) = 
            error("LsqExpansion2 only support Complex data currently")

function _exponential_expansion(f::Vector{<:Complex}, alg::LsqExpansion2; steps::Vector{Int})
    L = length(f)
    atol = alg.atol
    nitr = min(alg.n, L)
    xmax = length(f)
    xdata = collect(1:xmax)

    alps, lams = lsq_prony(f, 1)
    err = expansion_error(f, [alps; lams])
    best_n, best_err, best_alps, best_lams = 1, err, alps, lams

    for n in 1:nitr
        (alps, lams), err = _lsq_expansion_n(f, n, steps=steps)
        if err < best_err
            best_n, best_alps, best_lams, best_err = n, alps, lams, err
        end
        if best_err <= atol
            (alg.verbosity >= 1) && println("LsqExpansion2 converged in $n iterations, error is $best_err")
            return best_alps, best_lams
        end
        if n >= min(L-n+1, nitr)
            @warn "can not find a good approximation with L=$(L), n=$(alg.n), atol=$(atol), return with error $err"
            return best_alps, best_lams
        end
    end
    error("can not be here")
end

function _lsq_expansion_n(f::Vector{<:Complex}, p::Int; steps::Vector{Int})
    errs = []
    coeffs = []
    for step in steps
        f2 = f[1:step:end]
        (p > length(f2)÷2) && continue

        xs, lambdas = lsq_prony(f2, p)
        xs′ = @. xs * lambdas ^(1-1/step)
        lambdas′ = @. lambdas ^ (1/step)

        xs′, lambdas′, err = lsq_expansion_n(f, xs′, lambdas′)
        err = expansion_error(f, xs′, lambdas′)

        push!(errs, err)
        push!(coeffs, (xs′, lambdas′))
    end

    _, idx = findmin(identity, errs)
    return coeffs[idx], errs[idx]
end


function _lsq_predict_one(x::Integer, p::Vector{<:Real})
    n = div(length(p), 4)
    na = p[1:n] # norm
    pa = p[n+1:2*n] # phase
    nb = p[2*n+1:3*n] # norm
    pb = p[3*n+1:end] # phase

    T = eltype(p)
    rel = zero(T)
    img = zero(T)
    for j in 1:n
        tmp = na[j]*nb[j]^x
        phase = pa[j] + pb[j]*x
        rel += tmp * cos(phase)
        img += tmp * sin(phase)
    end
    return rel, img
end
function _lsq_predict(xmax::Integer, p::Vector{<:Real})
    @assert length(p) % 4 == 0
    res = [_lsq_predict_one(i, p) for i in 1:xmax]
    r = [[d[1] for d in res]; [d[2] for d in res]]
    return r
end
function _lsq_jacobian_one(x::Integer, p::Vector{<:Real})
    n = div(length(p), 4)
    na = p[1:n] # norm
    pa = p[n+1:2*n] # phase
    nb = p[2*n+1:3*n] # norm
    pb = p[3*n+1:end] # phase

    nab = @. na * nb^x
    pab = @. pa + pb*x
    nab_cos = @. nab * cos(pab)
    nab_sin = @. nab * sin(pab)
    # rel = sum(nab_cos)
    # img = sum(nab_sin)

    r_na = @. nab_cos / na
    r_nb = @. nab_cos / nb * x
    r_pa = -nab_sin
    r_pb = -nab_sin * x
    i_na = @. nab_sin / na
    i_nb = @. nab_sin / nb * x
    i_pa = nab_cos
    i_pb = nab_cos * x
    return [r_na; r_pa; r_nb; r_pb], [i_na; i_pa; i_nb; i_pb]
end
function _lsq_jacobian(xmax::Integer, p::Vector{<:Real})
    @assert length(p) % 4 == 0
    res = [_lsq_jacobian_one(i, p) for i in 1:xmax]
    r = [[d[1] for d in res]; [d[2] for d in res]]
    return permutedims(hcat(r...))
end
function lsq_expansion_n(f::Vector{<:Complex}, alphas::Vector{<:Complex}, lambdas::Vector{<:Complex})
    n = length(alphas)
    @assert n == length(lambdas)

    ydata = [real(f); imag(f)]
    xmax = length(f)
    xdata = collect(1:xmax)
    p0 = [norm.(alphas); angle.(alphas); norm.(lambdas); angle.(lambdas)]
    fit = curve_fit((x,p)->_lsq_predict(xmax,p), (x,p)->_lsq_jacobian(xmax,p), xdata, ydata, p0)

    p = fit.param
    alp = @. p[1:n] * exp(im*p[n+1:2*n])
    lam = @. p[2*n+1:3*n] * exp(im*p[3*n+1:4*n])
    
    err = expansion_error(f, [alp; lam])
    return alp, lam, err
end









function first_period(x::Vector{<:Real})
    idx = findfirst(i->!((x[i] > x[i-1]) ⊻ (x[i] > x[i+1])), 2:(length(x)-1))
    return isnothing(idx) ? length(x) : idx + 1
end

function cut(f::Vector, as, bs, alg::ExponentialExpansionAlgorithm2)
    p = sortperm(abs.(as), rev=true)
    as, bs = as[p], bs[p]
    N = length(as)

    (expansion_error(f, as, bs) >= alg.atol) && return as, bs
    while (expansion_error(f, as[1:N], bs[1:N]) < alg.atol) && (N > 1)
        N -= 1
    end
    return as[1:N+1], bs[1:N+1]
end


