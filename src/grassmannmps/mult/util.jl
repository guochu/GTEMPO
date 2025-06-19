function get_left_xy(left::MPSTensor, xj::MPSTensor, yj::MPSTensor)
    @tensor tmp1[1,5,4;2] := left[1,2,3] * yj[3,4,5]
    @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * xj[4,5,6]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end
    tmp3 = g_fuse(tmp2, 2)
    return tmp3
end
function get_xy_right(right::MPSTensor, xj::MPSTensor, yj::MPSTensor)
    @tensor tmp1[4; 1 2 5] := yj[1,2,3] * right[3,4,5]
    @tensor tmp2[4 1 2 5; 6] := xj[1,2,3] * tmp1[3,4,5,6]
    for (f1, f2) in fusiontrees(tmp2)
        coef = (isodd(f1.uncoupled[2].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end 
    tmp3 = g_fuse(tmp2, 3)
    return tmp3
end
