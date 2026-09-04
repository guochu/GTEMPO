#---------------------------------------------------------------
# @grassmann macro
#
# Same syntax as TensorOperations.@tensor, but contractions follow the
# fermionic (Grassmann) sign convention. This is achieved by running the
# standard TensorOperations parser with `GrassmannBackend()` inserted into
# every `tensoradd!`, `tensortrace!` and `tensorcontract!` call, which
# dispatches to the fermionic implementations for parity TensorMaps defined
# in tensoroperations.jl. No wrapping of the tensors is required: all input
# and output tensors are plain (Z2-graded) TensorMaps. Plain `@tensor`
# expressions are unaffected and keep the bosonic Z2Tensors semantics.
#---------------------------------------------------------------

"""
    @grassmann(tensor_expr)
    @grassmann [kwargs...] tensor_expr

Grassmann version of `TensorOperations.@tensor`: the allowed expressions and
keyword arguments (`order`, `opt`, `contractcheck`, `costcheck`, ...) are
exactly those of `@tensor`, but all tensors are plain Z2-graded `TensorMap`s
and the contractions, additions and traces are performed with the fermionic
(Grassmann) sign convention. This is equivalent to `@tensor` with
`backend = GrassmannBackend()`.

Examples:
```julia
@grassmann c[1,4,5;3] := a[1,2,3] * b[2,4,5]
@grassmann s = a[1,2] * b[2,1]           # scalar output
@grassmann c[1,2] += a[2,1]              # in-place addition
```
"""
macro grassmann(args...)
    isempty(args) && throw(ArgumentError("No arguments passed to `@grassmann`"))
    backend = Expr(:call, GlobalRef(GTEMPO, :GrassmannBackend))
    if length(args) == 1
        parser = TO.tensorparser(args[1], :backend => backend)
        return esc(parser(args[1]))
    end
    tensorexpr = args[end]
    kwargs = TO.parse_tensor_kwargs(args[1:(end - 1)])
    any(kw -> kw.first === :backend, kwargs) &&
        throw(ArgumentError("`@grassmann` does not support a custom `backend` keyword argument"))
    parser = TO.tensorparser(tensorexpr, kwargs..., :backend => backend)
    return esc(parser(tensorexpr))
end
