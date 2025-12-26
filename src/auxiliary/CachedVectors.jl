# v2

module CachedVectors

using Serialization
export CachedVector, all2disk!, change_sweeporder!, clear!, destory_copy!, @withlock

# suitable for DMRG-like sweep algorithms
# have a center point as reference to manage memory

# TODO:
# May can go further
# while !started[] is not necessary
# only need lock when read/write CachedVector.data 
# or maybe only when sub thread is changing CachedVector.data
# when sub thread is read/wirte file, no need lock 
# CANCLED: main thread may repeat ask sub thread to load/save


const CachedVectorPath = Ref{String}()
function __init__()
    current_pid = getpid()
    full_path = "CVTemporaryFiles/$current_pid/"
    
    # remove old folder with same PID
    if isdir(full_path)
        try
            rm(full_path; recursive=true, force=true)
        catch e
            @warn "Failed to clear old PID folder: $e"
        end
    end

    # create folder
    try
        mkpath(full_path)
        CachedVectorPath[] = full_path
    catch e
        error("Can't create folder: $e")
    end

    # clear folder when exit
    atexit() do
        if isassigned(CachedVectorPath) && !isempty(CachedVectorPath[])
            path_to_clean = CachedVectorPath[]
            if isdir(path_to_clean)
                try
                    rm(path_to_clean; recursive=true, force=true)
                    # println("Cleared the exiting folder -> $path_to_clean")
                catch e
                    @warn "Failed to clear folder $path_to_clean: $e"
                end
            end
        end
    end
end
get_cvpath() = CachedVectorPath[]


idx2file(idx::Int) = "$(get_cvpath())/cache_$idx"
function create_diskpath()
    paths = [p for p in readdir(get_cvpath(); join=true) if isdir(p)]
    if length(paths) == 0
        idx = 1
    else
        n = length(get_cvpath()) + 7
        idxs = sort(map(x->parse(Int, x[n:end]), paths))
        idx = findfirst(x->!insorted(x, idxs), 1:typemax(Int))
    end
    diskpath = idx2file(idx)
    mkpath(diskpath)
    return idx
end


mutable struct CachedVector{T} <: AbstractVector{T}
    data::Vector{Union{T,Missing}}
    locker::ReentrantLock
    diskpath::Int
    isleftsweep::Bool
    function CachedVector{T}(::UndefInitializer, L::Int, isleftsweep::Bool=true) where T
        diskpath = create_diskpath()
        data = Vector{T}(undef, L)
        obj = new{T}(data, ReentrantLock(), diskpath, isleftsweep)
        finalizer(obj) do c
            # @info "Object been GC-ed!"
            lock(c.locker) do
                c.diskpath != 0 && clear!(c)
            end
            c
        end
        return obj
    end
    function CachedVector(data::Vector{T}, isleftsweep::Bool=true) where T
        diskpath = create_diskpath()
        obj = new{T}(data, ReentrantLock(), diskpath, isleftsweep)
        finalizer(obj) do c
            # @info "Object been GC-ed!"
            lock(c.locker) do
                c.diskpath != 0 && clear!(c)
            end
            c
        end
        return obj
    end
end
Base.length(c::CachedVector) = length(c.data)
Base.firstindex(c::CachedVector) = 1
Base.lastindex(c::CachedVector) = length(c.data)
function diskpath(c::CachedVector)
    idx2file(c.diskpath)
end

macro withlock(obj, expr)
    return quote
        lock($(esc(obj)).locker) do
            $(esc(expr))
        end
    end
end

change_sweeporder!(c::CachedVector) = (c.isleftsweep = !(c.isleftsweep))
change_sweeporder!(c::CachedVector, isleftsweep::Bool) = (c.isleftsweep = isleftsweep)
center2idx(c::CachedVector, idx::Int) = c.isleftsweep ? (idx, idx+1) : (idx, idx-1)

function state(c::CachedVector{T}, i::Int) where T
    lock(c.locker) do
        if !isassigned(c.data, i)
            return :undef
        else
            if ismissing(c.data[i])
                return :indisk
            else
                return :inCPU
            end
        end
    end
end

function filename(c::CachedVector, idx::Int)
    joinpath(diskpath(c), "tmp$idx")
end
function Base.display(c::CachedVector{T}) where T
    println("CachedVector of $T")
    data = c.data
    for i in eachindex(data)
        println("$i ", state(c, i))
    end
end

# move data from/to disk
function _load_from_disk!(c::CachedVector{T}, idx::Int) where T
    @assert ismissing(c.data[idx])
    file = filename(c, idx)
    if isfile(file)
        dat = deserialize(file)
        c.data[idx] = dat
        rm(file)
        return dat
    else
        throw("index of $idx didn't in Disk")
    end
end

function _safe_serialize(path, obj)
    Serialization.serialize(path, obj)
    # get the correct permission: rw-r--r--
    try
        chmod(path, 0o644)
    catch e
        @warn "chmod failed" e
    end
    return nothing
end
function _save_to_disk!(c::CachedVector{T}, idx::Int) where T
    @assert c.data[idx] isa T
    _safe_serialize(filename(c, idx), c.data[idx])
    c.data[idx] = missing
end

# TODO: should be 1, but hleft/hright start from 2
function _auto_change_order!(c::CachedVector, idx::Int)
    if idx <= 2
        change_sweeporder!(c, true)
    elseif idx+2 > length(c)
        change_sweeporder!(c, false)
    end
end

function Base.getindex(c::CachedVector, idx::Int)
    local res
    locked_signal = Base.Event()
    lock(c.locker) do
        _auto_change_order!(c, idx)
        if isassigned(c.data, idx)
            res = if ismissing(c.data[idx])
                _load_from_disk!(c, idx)
            else
                c.data[idx]
            end
        else
            error("UndefRefError: access to undefined reference")
        end

        Threads.@spawn begin
            lock(c.locker) do
                notify(locked_signal)
                all2disk!(c; except = center2idx(c, idx))
            end
        end
    end
    wait(locked_signal)
    return res
end

function Base.setindex!(c::CachedVector{T}, data::T, idx::Int) where T
    locked_signal = Base.Event()
    lock(c.locker) do
        _auto_change_order!(c, idx)
        # println("$(Threads.threadid()) Set $idx")
        if isassigned(c.data, idx) && ismissing(c.data[idx])
            rm(filename(c, idx))
        end
        c.data[idx] = data

        s = Threads.@spawn begin
            lock(c.locker) do
                notify(locked_signal)
                all2disk!(c; except = center2idx(c, idx))
            end
        end    
    end
    wait(locked_signal)

    # fetch(s)
    # println("Hello")
end

function all2disk!(c::CachedVector{T}; except::NTuple{2,Int}) where T
    lock(c.locker) do
        # println("$(Threads.threadid()) all2disk 2")

        for i in eachindex(c.data)
            if isassigned(c.data, i)
                if i in except
                    if ismissing(c.data[i])
                        _load_from_disk!(c, i)
                    end
                elseif !ismissing(c.data[i])
                    _save_to_disk!(c, i)
                end
            end
        end
    end
end

function clear!(c::CachedVector)
    !ispath(get_cvpath()) && return
    lock(c.locker) do
        for i in eachindex(c.data)
            if isassigned(c.data, i)
                if ismissing(c.data[i])
                    rm(filename(c, i))
                else
                    c.data[i] = missing
                end
            end
        end
        rm(diskpath(c))
        c.diskpath = 0
    end
end

function destory_copy!(x::CachedVector{T}, y::CachedVector{T}) where T
    @assert length(x) == length(y)
    clear!(x)
    lock(x.locker) do
        lock(y.locker) do
            x.data .= y.data
            x.diskpath = y.diskpath
            x.isleftsweep = y.isleftsweep
            y.data .= missing
            y.diskpath = 0
        end
    end
    return x
end
function destory_copy!(x::Vector{T}, y::CachedVector{T}) where T
    @assert length(x) == length(y)
    for idx in 1:length(x)
        if isassigned(y.data, idx)
            x[idx] = if ismissing(y.data[idx])
                _load_from_disk!(y, idx)
            else
                y.data[idx]
            end
        else
            error("UndefRefError: access to undefined reference")
        end
    end
    clear!(y)
    return x
end

end

# using .CachedVectors
@reexport using .CachedVectors

CachedVectors.destory_copy!(x::Vector{T}, y::Vector{T}) where T = copy!(x, y)

const DefaultUseCache = false
# maybe it's more convenient to decide whether use cache according to the memory needed
const MemoryThreshold = 2^30 * 16 # 16 GiB


function test_CachedVector()
    diskpath = CachedVectors.diskpath
    set(c, i) = c[i] = fill(i, 100, 500)
    check(c, i) = all(c[i] .== i)

    c1 = CachedVector{Matrix}(undef, 10)
    set(c1, 1)
    @withlock c1 @assert [p for p in readdir(diskpath(c1))] == []
    set(c1, 8)
    @withlock c1 @assert [p for p in readdir(diskpath(c1))] == ["tmp1"]
    set(c1, 2)
    @withlock c1 @assert sort([p for p in readdir(diskpath(c1))]) == ["tmp1", "tmp8"]
    set(c1, 1)
    @withlock c1 @assert sort([p for p in readdir(diskpath(c1))]) == ["tmp8"]

    @assert check(c1, 1)
    @assert check(c1, 2)
    @assert check(c1, 8)

    all2disk!(c1, except = (2,7))
    @withlock c1 @assert sort([p for p in readdir(diskpath(c1))]) == ["tmp1", "tmp8"]

    change_sweeporder!(c1)
    set(c1, 3)
    @withlock c1 @assert [p for p in readdir(diskpath(c1))] == ["tmp1", "tmp8"]
    change_sweeporder!(c1, false)
    set(c1, 7)
    @withlock c1 @assert [p for p in readdir(diskpath(c1))] == ["tmp1", "tmp2", "tmp3", "tmp8"]
    set(c1, 2)
    @withlock c1 @assert sort([p for p in readdir(diskpath(c1))]) == ["tmp1", "tmp7", "tmp8"]

    change_sweeporder!(c1, true)
    set(c1, 7)
    @withlock c1 @assert sort([p for p in readdir(diskpath(c1))]) == ["tmp1", "tmp2", "tmp3"]


    change_sweeporder!(c1)

    erridx = []
    for i in 1:10
        try
            c1[i]
        catch e
            push!(erridx, i)
        end
    end
    @assert erridx == [4,5,6,9,10]

    for i in 1:10
        set(c1, i)
    end
    for i in 1:10
        @assert check(c1, i)
    end

    @withlock c1 @assert sort([p for p in readdir(diskpath(c1))]) == ["tmp$i" for i in 1:8]

    c2 = CachedVector{Matrix}(undef, 10)
    set(c2, 1)
    set(c2, 8)
    set(c2, 2)
    set(c2, 1)
    @withlock c2 @assert sort([p for p in readdir(diskpath(c2))]) == ["tmp8"]

    destory_copy!(c2, c1)
    v = Vector{Matrix}(undef, 10)
    destory_copy!(v, c2)
    for i in 1:10
        @assert check(v, i)
    end


    # clear!(c1)
end
# test_CachedVector()


