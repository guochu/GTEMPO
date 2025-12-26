
function parint_mult(xs::GrassmannMPS...; cidx::Vector{Int}, trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
    Lxs = length(xs[1])
    (unique(length.(xs)) == [Lxs,]) || throw(DimensionMismatch())
    check_contract_idx(Lxs, cidx)

    Lz = Lxs-2*length(cidx)
    A = xs[1][1]
	B = bondtensortype(spacetype(A), real(scalartype(A)))
	z = GrassmannMPS(similar(xs[1].data, Lz), Vector{Union{Missing, B}}(missing, Lz), Ref(1.0))
    left = isomorphism(scalartype(xs[1]), fuse(space_l.(xs)...), ⊗(space_l.(xs)...) )
    right = isomorphism(scalartype(xs[1]), ⊗(space_r.(xs)...)', fuse(space_r.(xs)...))

    ixs = 1
    iz = 1
	while true
		rt = @elapsed if insorted(ixs, cidx)
            left = left_m(left, GrassmannTransferMatrix((ixs+1)÷2, xs...))
            ixs += 2
		else
            tmp = get_left_below(left, getindex.(xs, ixs)...)
            z[iz], left = leftorth!(tmp, alg=QR())
            iz += 1
            ixs += 1
		end
        _renormalize!(z, left, false)
        (verbosity >= 2) && println("$ixs / $Lxs cost $rt Seconds, with bond tensor of dimension $(dim(space(left, 1)))")
        (ixs > Lxs) && break
	end
    @tensor tmp[1 2;4] := z[end][1 2 3] * left_right(left, right)[3 4]
    z[end] = tmp

    # z.svectors[1] = Diagonal(id(space_l(z[1])))
	# z.svectors[end] = Diagonal(id(space_r(z[end])'))
	z.svectors[1] = DiagonalTensorMap{Float64}(ones, space_l(z[1]) )
	z.svectors[end] = DiagonalTensorMap{Float64}(ones, space_r(z[end])' )

    setscaling!(z, *(scaling.(xs)...) ^ (length(xs[1]) / length(z)) * scaling(z))
    (verbosity >= 2) && println("bond dimension of intermediate GMPS: ", bond_dimension(z))
    _rightorth!(z, SVD(), trunc, false, verbosity)
    return z
end


# multiply n GMPS together
# left to right
function my_mult(xs::GrassmannMPS...; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
    Lxs = length(xs[1])
    (unique(length.(xs)) == [Lxs,]) || throw(DimensionMismatch("unique($(length.(xs))) != [$Lxs,]"))
    res = copy(xs[1])

    left = isomorphism(scalartype(xs[1]), fuse(space_l.(xs)...), ⊗(space_l.(xs)...) )
    right = isomorphism(scalartype(xs[1]), ⊗(space_r.(xs)...)', fuse(space_r.(xs)...))

    tmp4 = get_left_below(left, getindex.(xs, 1)...)
    for i in 1:Lxs-1
        q, r = leftorth!(tmp4, alg = QR())
        res[i] = q
        _renormalize!(res, r, false)
        tmp4 = get_left_below(r, getindex.(xs, i+1)...)
    end
    tmp = left_below__right(tmp4, right)
    res[end] = tmp
    _rightorth!(res, SVD(), trunc, false, verbosity)
    setscaling!(res, *(scaling(res), scaling.(xs[2:end])...))
    return res
end
# right to left
function my_mult2(xs::GrassmannMPS...; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
    nx = length(xs)
    Lxs = length(xs[1])
    (unique(length.(xs)) == [Lxs,]) || throw(DimensionMismatch("unique($(length.(xs))) != [$Lxs,]"))
    res = copy(xs[1])

    left = isomorphism(scalartype(xs[1]), fuse(space_l.(xs)...), ⊗(space_l.(xs)...) )
    right = isomorphism(scalartype(xs[1]), ⊗(space_r.(xs)...)', fuse(space_r.(xs)...))

    tmp4 = get_below_right(right, getindex.(xs, Lxs)...)
    for i in Lxs:-1:2
        l, q = rightorth(tmp4, alg=LQ())
        res[i] = permute(q, (1,2), (3,))

        _renormalize!(res, l, false)
        tmp4 = get_below_right(l, getindex.(xs, i-1)...)
    end
    tmp = left__below_right(left, tmp4)
    res[1] = tmp
    _leftorth!(res, SVD(), trunc, false, verbosity)
    setscaling!(res, *(scaling(res), scaling.(xs[2:end])...))
    return res
end

