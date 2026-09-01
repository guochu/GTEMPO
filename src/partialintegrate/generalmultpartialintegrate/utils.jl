
# left (z,), (xs...)
# left_below (z, phy), (xs...)
# right (xs...), (z,)
# below_right (xs...), (phy, z)

function get_left_below(left::AbstractParityTensorMap{<:Number, 1, 2}, xj::MPSTensor, yj::MPSTensor)
    @tensor tmp1[1 5 4;2] := left[1 2 3] * yj[3 4 5]
    @tensor tmp2[1 4 6;7 5] := tmp1[1 5 4 2] * xj[2 6 7]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end
    tmp3 = g_fuse(tmp2, 2)
    return tmp3
end
function get_left_below(left::AbstractParityTensorMap{<:Number, 1, 3}, xj::MPSTensor, yj::MPSTensor, zj::MPSTensor)
    @tensor tmp1[1 2 5 6;3] := left[1 2 3 4] * zj[4 5 6]
    @tensor tmp2[1 5 7 8 6;2] := tmp1[1 2 5 6 3] * yj[3 7 8]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f1.uncoupled[2].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end
    tmp3 = g_fuse(tmp2, 2)

    @tensor tmp4[1 57 9;a 8 6] := tmp3[1 57 8 6 2] * xj[2 9 a]
    for (f1, f2) in fusiontrees(tmp4)
        coef = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp4[f1, f2])
        end
    end
    tmp5 = g_fuse(tmp4, 2)
    return tmp5
end
function get_left_below(left::AbstractParityTensorMap{<:Number, 1, 4}, xj::MPSTensor, yj::MPSTensor, zj::MPSTensor, mj::MPSTensor)
    @tensor tmp1[1 2 3 6 7;4] := left[1 2 3 4 5] * mj[5 6 7]
    @tensor tmp2[1 2 8 6 9 7;3] := tmp1[1 2 3 6 7 4] * zj[4 8 9]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end
    tmp3 = g_fuse(tmp2, 3)

    @tensor tmp4[1 a 86 b 9 7;2] := tmp3[1 2 86 9 7 3] * yj[3 a b]
    for (f1, f2) in fusiontrees(tmp4)
        coef = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp4[f1, f2])
        end
    end
    tmp5 = g_fuse(tmp4, 2)

    @tensor tmp6[1 c a86;d b 9 7] := tmp5[1 a86 b 9 7 2] * xj[2 c d]
    for (f1, f2) in fusiontrees(tmp6)
        coef = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp6[f1, f2])
        end
    end
    tmp7 = g_fuse(tmp6, 2)
    return tmp7
end
function get_left_below(left::AbstractParityTensorMap{<:Number, 1, 5}, xj::MPSTensor, yj::MPSTensor, zj::MPSTensor, mj::MPSTensor, nj::MPSTensor)
    @tensor tmp1[1 2 3 4 7 8;5] := left[1 2 3 4 5 6] * nj[6 7 8]
    @tensor tmp2[1 2 3 9 7 a 8;4] := tmp1[1 2 3 4 7 8 5] * mj[5 9 a]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end
    tmp3 = g_fuse(tmp2, 4)

    @tensor tmp4[1 2 b 97 c a 8;3] := tmp3[1 2 3 97 a 8 4] * zj[4 b c]
    for (f1, f2) in fusiontrees(tmp4)
        coef = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp4[f1, f2])
        end
    end
    tmp5 = g_fuse(tmp4, 3)

    @tensor tmp6[1 d b97 e c a 8;2] := tmp5[1 2 b97 c a 8 3] * yj[3 d e]
    for (f1, f2) in fusiontrees(tmp6)
        coef = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp6[f1, f2])
        end
    end
    tmp7 = g_fuse(tmp6, 2)

    @tensor tmp8[1 f db97;g e c a 8] := tmp7[1 db97 e c a 8 2] * xj[2 f g]
    for (f1, f2) in fusiontrees(tmp8)
        coef = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp8[f1, f2])
        end
    end
    tmp9 = g_fuse(tmp8, 2)
    return tmp9
end
function get_left_below(left::AbstractParityTensorMap{<:Number, 1, 6}, xj::MPSTensor, yj::MPSTensor, zj::MPSTensor, mj::MPSTensor, nj::MPSTensor, pj::MPSTensor)
    @tensor tmp1[1 2 3 4 5 8 9;6] := left[1 2 3 4 5 6 7] * pj[7 8 9]
    @tensor tmp2[1 2 3 4 a 8 b 9;5] := tmp1[1 2 3 4 5 8 9 6] * nj[6 a b]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end
    tmp3 = g_fuse(tmp2, 5)

    @tensor tmp4[1 2 3 c a8 d b 9;4] := tmp3[1 2 3 4 a8 b 9 5] * mj[5 c d]
    for (f1, f2) in fusiontrees(tmp4)
        coef = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp4[f1, f2])
        end
    end
    tmp5 = g_fuse(tmp4, 4)

    @tensor tmp6[1 2 e ca8 f d b 9;3] := tmp5[1 2 3 ca8 d b 9 4] * zj[4 e f]
    for (f1, f2) in fusiontrees(tmp6)
        coef = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp6[f1, f2])
        end
    end
    tmp7 = g_fuse(tmp6, 3)

    @tensor tmp8[1 g eca8 h f d b 9;2] := tmp7[1 2 eca8 f d b 9 3] * yj[3 g h]
    for (f1, f2) in fusiontrees(tmp8)
        coef = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp8[f1, f2])
        end
    end
    tmp9 = g_fuse(tmp8, 2)

    @tensor tmp10[1 i geca8;j h f d b 9] := tmp9[1 geca8 h f d b 9 2] * xj[2 i j]
    for (f1, f2) in fusiontrees(tmp10)
        coef = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp10[f1, f2])
        end
    end
    tmp11 = g_fuse(tmp10, 2)
    return tmp11
end
# can be implemented for more xs


function get_below_right(right::AbstractParityTensorMap{<:Number, 2, 1}, xj::MPSTensor, yj::MPSTensor)
    @tensor tmp1[2; 4 5 3] := yj[4 5 1] * right[1 2 3]
    @tensor tmp2[4 6; 7 5 3] := xj[6 7 2] * tmp1[2 4 5 3]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end 
    tmp3 = g_fuse(tmp2, 3)
    return tmp3
end
function get_below_right(right::AbstractParityTensorMap{<:Number, 3, 1}, xj::MPSTensor, yj::MPSTensor, zj::MPSTensor)
    @tensor tmp1[2; 5 6 3 4] := zj[5 6 1] * right[1 2 3 4]
    @tensor tmp2[3;5 7 6 8 4] := yj[7 8 2] * tmp1[2 5 6 3 4]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end 
    tmp3 = g_fuse(tmp2, 4)

    @tensor tmp4[5 7 9; 68 a 4] := xj[9 a 3] * tmp3[3 5 7 68 4]
    for (f1, f2) in fusiontrees(tmp4)
        coef = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp4[f1, f2])
        end
    end 
    tmp5 = g_fuse(tmp4, 4)
    return tmp5
end
function get_below_right(right::AbstractParityTensorMap{<:Number, 4, 1}, xj::MPSTensor, yj::MPSTensor, zj::MPSTensor, mj::MPSTensor)
    @tensor tmp1[2; 6 7 3 4 5] := mj[6 7 1] * right[1 2 3 4 5]
    @tensor tmp2[3;6 8 7 9 4 5] := zj[8 9 2] * tmp1[2 6 7 3 4 5]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end 
    tmp3 = g_fuse(tmp2, 4)

    @tensor tmp4[4; 6 8 a 79 b 5] := yj[a b 3] * tmp3[3 6 8 79 4 5]
    for (f1, f2) in fusiontrees(tmp4)
        coef = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[4].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp4[f1, f2])
        end
    end 
    tmp5 = g_fuse(tmp4, 5)

    @tensor tmp6[6 8 a c; b79 d 5] := xj[c d 4] * tmp5[4 6 8 a b79 5]
    for (f1, f2) in fusiontrees(tmp6)
        coef = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp6[f1, f2])
        end
    end 
    tmp7 = g_fuse(tmp6, 5)
    return tmp7
end
function get_below_right(right::AbstractParityTensorMap{<:Number, 5, 1}, xj::MPSTensor, yj::MPSTensor, zj::MPSTensor, mj::MPSTensor, nj::MPSTensor)
    @tensor tmp1[2;7 8 3 4 5 6] := nj[7 8 1] * right[1 2 3 4 5 6]
    @tensor tmp2[3;7 9 8 a 4 5 6] := mj[9 a 2] * tmp1[2 7 8 3 4 5 6]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end 
    tmp3 = g_fuse(tmp2, 4)

    @tensor tmp4[4;7 9 b a8 c 5 6] := zj[b c 3] * tmp3[3 7 9 a8 4 5 6]
    for (f1, f2) in fusiontrees(tmp4)
        coef = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[4].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp4[f1, f2])
        end
    end 
    tmp5 = g_fuse(tmp4, 5)

    @tensor tmp6[5;7 9 b d ca8 e 6] := yj[d e 4] * tmp5[4 7 9 b ca8 5 6]
    for (f1, f2) in fusiontrees(tmp6)
        coef = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[5].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp6[f1, f2])
        end
    end 
    tmp7 = g_fuse(tmp6, 6)

    @tensor tmp8[7 9 b d f;eca8 g 6] := xj[f g 5] * tmp7[5 7 9 b d eca8 6]
    for (f1, f2) in fusiontrees(tmp8)
        coef = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp8[f1, f2])
        end
    end 
    tmp9 = g_fuse(tmp8, 6)
    return tmp9
end
function get_below_right(right::AbstractParityTensorMap{<:Number, 6, 1}, xj::MPSTensor, yj::MPSTensor, zj::MPSTensor, mj::MPSTensor, nj::MPSTensor, pj::MPSTensor)
    @tensor tmp1[2;8 9 3 4 5 6 7] := pj[8 9 1] * right[1 2 3 4 5 6 7]
    @tensor tmp2[3;8 a 9 b 4 5 6 7] := nj[a b 2] * tmp1[2 8 9 3 4 5 6 7]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end 
    tmp3 = g_fuse(tmp2, 4)

    @tensor tmp4[4;8 a c b9 d 5 6 7] := mj[c d 3] * tmp3[3 8 a b9 4 5 6 7]
    for (f1, f2) in fusiontrees(tmp4)
        coef = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[4].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp4[f1, f2])
        end
    end 
    tmp5 = g_fuse(tmp4, 5)

    @tensor tmp6[5;8 a c e db9 f 6 7] := zj[e f 4] * tmp5[4 8 a c db9 5 6 7]
    for (f1, f2) in fusiontrees(tmp6)
        coef = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[5].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp6[f1, f2])
        end
    end 
    tmp7 = g_fuse(tmp6, 6)

    @tensor tmp8[6;8 a c e g fdb9 h 7] := yj[g h 5] * tmp7[5 8 a c e fdb9 6 7]
    for (f1, f2) in fusiontrees(tmp8)
        coef = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[6].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp8[f1, f2])
        end
    end 
    tmp9 = g_fuse(tmp8, 7)

    @tensor tmp10[8 a c e g i;hfdb9 j 7] := xj[i j 6] * tmp9[6 8 a c e g hfdb9 7]
    for (f1, f2) in fusiontrees(tmp10)
        coef = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp10[f1, f2])
        end
    end 
    tmp11 = g_fuse(tmp10, 7)
    return tmp11
end
# can be implemented for more xs




# return a MPSTensor
function left__below_right(left::AbstractParityTensorMap{<:Number, 1, N}, below_right::AbstractParityTensorMap{<:Number, N, 2}) where N
    TO.tensorcontract(left, ((1,), ntuple(i->i+1, N)), false, 
                    below_right, (ntuple(i->N+1-i, N), (N+1,N+2)), false, ((1,2), (3,)))
end
function left_below__right(left_below::AbstractParityTensorMap{<:Number, 2, N}, right::AbstractParityTensorMap{<:Number, N, 1}) where N
    TO.tensorcontract(left_below, ((1,2), ntuple(i->i+2, N)), false, 
                    right, (ntuple(i->N+1-i, N), (N+1,)), false, ((1,2), (3,)))
end

# return a bond Tensor
function left_right(left::AbstractParityTensorMap{<:Number, 1, N}, right::AbstractParityTensorMap{<:Number, N, 1}) where N
    TO.tensorcontract(left, ((1,), ntuple(i->i+1, N)), false, 
                    right, (ntuple(i->N+1-i, N), (N+1,)), false, ((1,), (2,)))
end
# return a scalar Tensor
function leftright(left::AbstractParityTensorMap{<:Number, 1, N}, right::AbstractParityTensorMap{<:Number, N, 1}) where N
    TO.tensorcontract(left, ((), ntuple(i->i, N+1)), false, 
                    right, (ntuple(i->N+2-i, N+1), ()), false, ((), ()))
end

# return a left boundary Tensor
function left_below_above(left_below::AbstractParityTensorMap{<:Number, 2, N}, above::MPSTensor) where N
    TO.tensorcontract(left_below, (ntuple(i->i+2, N), (1,2)), false, 
                    above, ((1,2), (3,)), true, ((N+1,), ntuple(i->i, N)))
end
# return a right boundary Tensor
function above_below_right(above::MPSTensor, below_right::AbstractParityTensorMap{<:Number, N, 2}) where N
    TO.tensorcontract(above, ((1,), (2,3)), true, 
                    below_right, ((N+1,N+2), ntuple(i->i, N)), false, (ntuple(i->i+1,N), (1,)))
end


# need to remove the scaling in transfer matrix
function rescaling(m::GrassmannTransferMatrix)
    states = map(x -> getproperty.(x, :data), m.states)
    GrassmannTransferMatrix(states)
end
# need to reverse the column index, when convert TensorMap to GrassmannTensorMap
function left_m(left::AbstractParityTensorMap{<:Number, 1, N}, m::GrassmannTransferMatrix) where N
    left = deepcopy(left)
    p = (1,ntuple(i->N+2-i,N)...)
    swaps = TK.permutation2swaps(p)
    for (f1, f2) in fusiontrees(left)
        ns = getproperty.([f1.uncoupled..., dual.(f2.uncoupled)...], :n)
        coef = 1
    	for s in swaps
    		coef = (isodd(ns[s]) && isodd(ns[s+1])) ? -coef : coef
            ns[s], ns[s+1] = ns[s+1], ns[s]
    	end
        if coef != 1
            lmul!(coef, left[f1, f2])
        end
    end
    left = (GrassmannTensorMap(left) * rescaling(m)).data
    for (f1, f2) in fusiontrees(left)
        ns = getproperty.([f1.uncoupled..., dual.(f2.uncoupled)...], :n)
        coef = 1
    	for s in swaps
    		coef = (isodd(ns[s]) && isodd(ns[s+1])) ? -coef : coef
            ns[s], ns[s+1] = ns[s+1], ns[s]
    	end
        if coef != 1
            lmul!(coef, left[f1, f2])
        end
    end
    return left

    # another implementation
    # left = permute(left, (ntuple(i->i,N), (N+1,)))
    # left = permute(GrassmannTensorMap(left), ((1,),ntuple(i->i+1,N)))
    # left = left * m
    # left = permute(left, (ntuple(i->i,N), (N+1,)))
    # permute(left.data, ((1,),ntuple(i->i+1,N)))
    # return left
end
function m_right(m::GrassmannTransferMatrix, right::AbstractParityTensorMap{<:Number, N, 1}) where N
    (rescaling(m) * GrassmannTensorMap(right)).data
end


