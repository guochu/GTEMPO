TK.permute(m::AbstractArray, perm) = PermutedDimsArray(m, perm)
TK.permute(m::AbstractArray, left, right) = permute(m, (left..., right...))

function _eye(::Type{T}, m::Int, n::Int) where {T<:Number}
  r = zeros(T, m, n)
  for i in 1:min(m, n)
      r[i, i] = 1
  end
  return r
end
_eye(::Type{T}, d::Int) where {T<:Number} = _eye(T, d, d)
_eye(d::Int) = _eye(Float64, d)

# diag(m::AbstractVector{<: Number}) = LinearAlgebra.diagm(0=>m)

function _group_extent(extent::NTuple{N, Int}, idx::NTuple{N1, Int}) where {N, N1}
    ext = Vector{Int}(undef, N1)
    l = 0
    for i=1:N1
        ext[i] = prod(extent[(l+1):(l+idx[i])])
        l += idx[i]
    end
    return NTuple{N1, Int}(ext)
end


function tie(a::AbstractArray{T, N}, axs::NTuple{N1, Int}) where {T, N, N1}
    (sum(axs) != N) && error("total number of axes should equal to tensor rank.")
    return reshape(a, _group_extent(size(a), axs))
end

function Base.kron(a::AbstractArray{Ta, N}, b::AbstractArray{Tb, N}) where {Ta<:Number, Tb<:Number, N}
    N == 0 && error("empty tensors.")
    sa = size(a)
    sb = size(b)
    sc = Tuple(sa[i]*sb[i] for i=1:N)
    c = Array{promote_type(Ta, Tb), N}(undef, sc)
    ranges = Vector{UnitRange{Int}}(undef, N)
    for index in CartesianIndices(a)
        # ranges[1] = (index[1]*sb[1]+1):(index[1]+1)*sb[1]
        for j = 1:N
            ranges[j] = ((index[j]-1)*sb[j]+1):(index[j]*sb[j])
        end
        c[ranges...] = a[index]*b
    end
    return c
end

function stable_svd!(a::StridedArray{T, 2}, workspace::AbstractVector{T}) where T
    if length(workspace) < length(a)
        resize!(workspace, length(a))
    end
    ac = reshape(view(workspace, 1:length(a)), size(a))
    copy!(ac, a)
    try
        return TK.MatrixAlgebra.svd!(ac, TK.SDD())
    catch
        return TK.MatrixAlgebra.svd!(a, SVD())
    end
end

function TK.tsvd!(a::StridedArray{T, 2}, workspace::AbstractVector{T}=similar(a, length(a)); trunc::TruncationScheme=NoTruncation()) where {T}
    u, s, v = stable_svd!(a, workspace)
    d_old = length(s)
    s, err = _truncate!(s, trunc)
    d = length(s)
    if d == d_old
        return u, s, v, err
    else
        return u[:, 1:d], s, v[1:d, :], err
    end
end

function TK.tsvd!(a::StridedArray{T, N}, left::NTuple{N1, Int}, right::NTuple{N2, Int}, workspace::AbstractVector{T}=similar(a, length(a)); 
    trunc::TruncationScheme=NoTruncation()) where {T <: Number, N, N1, N2}
    (N == N1 + N2) || throw(DimensionMismatch())
    if length(workspace) <= length(a)
        resize!(workspace, length(a))
    end
    dim = (left..., right...)
    b = permute(a, dim)
    shape_b = size(b)
    ushape = shape_b[1:N1]
    vshape = shape_b[(N1+1):end]
    s1 = prod(ushape)
    s2 = prod(vshape)
    # u, s, v = F.U, F.S, F.Vt
    bmat = copy!(reshape(view(workspace, 1:length(a)), s1, s2), reshape(b, s1, s2))
    u, s, v, err = tsvd!(bmat, reshape(a, length(a)), trunc=trunc)
    md = length(s)
    return reshape(u, (ushape..., md)), s, reshape(v, (md, vshape...)), err
end

function tqr!(a::StridedMatrix) 
    q, r = TK.MatrixAlgebra.leftorth!(a, QR(), 0)
    return q, r
end

function tlq!(a::StridedMatrix) 
    l, q = TK.MatrixAlgebra.rightorth!(a, LQ(), 0)
    return l, q
end

"""
    qr(a::AbstractArray{T, N}, axs::Tuple{NTuple{N1, Int}, NTuple{N2, Int}}) where {T, N, N1, N2}
QR decomposition of QTensor a, by joining axs to be the second dimension
"""
function tqr!(a::AbstractArray{T, N}, left::NTuple{N1, Int}, right::NTuple{N2, Int}, workspace::AbstractVector{T}=similar(a, length(a))) where {T<:Number, N, N1, N2}
    (N == N1 + N2) || throw(DimensionMismatch())
    if length(workspace) <= length(a)
        resize!(workspace, length(a))
    end
    newindex = (left..., right...)
    a1 = permute(a, newindex)
    shape_a = size(a1)
    dimu = shape_a[1:N1]
    s1 = prod(dimu)
    dimv = shape_a[(N1+1):end]
    s2 = prod(dimv)
    bmat = copyto!(reshape(view(workspace, 1:length(a)), s1, s2), reshape(a1, s1, s2))
    # F = LinearAlgebra.qr!(bmat)
    # u = Base.Matrix(F.Q)
    # v = Base.Matrix(F.R)
    u, v = tqr!(bmat)
    s = size(v, 1)
    return reshape(u, dimu..., s), reshape(v, s, dimv...)
end

function tlq!(a::AbstractArray{T, N}, left::NTuple{N1, Int}, right::NTuple{N2, Int}, workspace::AbstractVector{T}=similar(a, length(a))) where {T<:Number, N, N1, N2}
    (N == N1 + N2) || throw(DimensionMismatch())
    if length(workspace) <= length(a)
        resize!(workspace, length(a))
    end
    newindex = (left..., right...)
    a1 = permute(a, newindex)
    shape_a = size(a1)
    dimu = shape_a[1:N1]
    s1 = prod(dimu)
    dimv = shape_a[(N1+1):end]
    s2 = prod(dimv)
    bmat = copyto!(reshape(view(workspace, 1:length(a)), s1, s2), reshape(a1, s1, s2))
    # F = LinearAlgebra.lq!(bmat)
    # u = Matrix(F.L)
    # v = Matrix(F.Q)
    u, v = tlq!(bmat)
    s = size(v, 1)
    return reshape(u, dimu..., s), reshape(v, s, dimv...)
end


# truncation
_truncate!(v::AbstractVector{<:Real}, trunc::NoTruncation, p::Real=2) = v, 0.
function _truncate!(v::AbstractVector{<:Real}, trunc::DMRG.TruncationDimCutoff, p::Real=2)
    sca = norm(v, p)
    dtrunc = findlast(Base.Fix2(>, sca * trunc.ϵ), v)
    if isnothing(dtrunc)
        dtrunc = trunc.add_back
    end
    dtrunc = min(trunc.D, dtrunc)
    err = norm(v[dtrunc+1:end], p)
    return  v[1:dtrunc], err / sca
end
