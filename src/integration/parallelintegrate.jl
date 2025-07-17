using Base.Threads
### calculate greens' functions

get_size_1(total_itr::Int, n_threads::Int, thread_id::Int) = div(total_itr * thread_id, n_threads)
get_size_2(total_itr::Int, n_threads::Int, thread_id::Int) = div(total_itr * (thread_id + 1), n_threads)

"""
    The inner function f is assumed to be Zero Based
"""
function parallel_run(total_itr::Int, n_threads::Int, f::Function, args...; kwargs...)
    if (n_threads > 1)
        # nave, resi, jobs = partition_jobs(total_itr, n_threads)
        Threads.@threads for thread_id in 0:(n_threads-1)
            ist = get_size_1(total_itr, n_threads, thread_id)+1
            ifn = get_size_2(total_itr, n_threads, thread_id) 
            f(ist, ifn, args...; kwargs...)
        end
    else
        f(1, total_itr, args...; kwargs...)
    end
end


function parallel_integrate(lattice::AbstractGrassmannLattice, x0::Vector{<:GrassmannMPS}, x1::GrassmannMPS...; kwargs...)
    Zs = zeros(scalartype(x0[1]), length(x0))
    function _f(ist::Int, ifn::Int, r, lat, A0, args...; kws...)
        for i in ist:ifn
            r[i] = integrate(lat, A0[i], args...; kws...)
        end
    end
    parallel_run(length(x0), Threads.nthreads(), _f, Zs, lattice, x0, x1...; kwargs...)
    return sum(Zs)
end
parallel_integrate(lattice::AbstractGrassmannLattice, x0::GrassmannMPS, x1::GrassmannMPS...; kwargs...) = integrate(lattice, x0, x1...; kwargs...)