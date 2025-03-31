
function integrate(x::FockMPS)
	L = length(x)
	sca = scaling(x)
	v = dropdims(sum(x[L], dims=2), dims=(2,3)) * sca
	for i in L-1:-1:1
		tmp = dropdims(sum(x[i], dims=2), dims=2) * sca
		v = tmp * v
	end
	return only(v)
end