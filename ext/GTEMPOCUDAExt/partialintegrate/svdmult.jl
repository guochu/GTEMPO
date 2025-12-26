
function cu_parint_mult(xs::GrassmannMPS...; cidx::Vector{Int}, trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
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
    left = tocu(left)
	while true
		rt = @elapsed if insorted(ixs, cidx)
            left = left_m(left, tocu(GrassmannTransferMatrix((ixs+1)÷2, xs...)))
            ixs += 2
		else
            tmp = get_left_below(left, tocu.(getindex.(xs, ixs))...)
            zi, left = leftorth!(tmp, alg=QR())
            z[iz] = fromcu(zi)
            iz += 1
            ixs += 1
		end
        _renormalize!(z, left, false)
        (verbosity >= 2) && println("$ixs / $Lxs cost $rt Seconds, with bond tensor of dimension $(dim(space(left, 1)))")
        (ixs > Lxs) && break
	end
    left = fromcu(left)
    @tensor tmp[1 2;4] := z[end][1 2 3] * left_right(left, right)[3 4]
    z[end] = tmp

    # z.svectors[1] = Diagonal(id(space_l(z[1])))
	# z.svectors[end] = Diagonal(id(space_r(z[end])'))
    z.svectors[1] = DiagonalTensorMap{Float64}(ones, space_l(z[1]) )
	z.svectors[end] = DiagonalTensorMap{Float64}(ones, space_r(z[end])' )

    setscaling!(z, *(scaling.(xs)...) ^ (length(xs[1]) / length(z)) * scaling(z))
    (verbosity >= 2) && println("bond dimension of intermediate GMPS: ", bond_dimension(z))
    _cu_rightorth!(z, SVD(), trunc, false, verbosity)
    return z
end

