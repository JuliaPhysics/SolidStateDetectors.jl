# Internal functionality of SolidStateDetectors, not part of the public API:

const _CPUArrayLike = Union{Array, SubArray{<:Any,<:Any,<:Array}}
const _LinIndexCPUArrayLike = Union{Array, SubArray{<:Any,<:Any,<:Array,<:Any,true}}
const _NonLinIndexView = SubArray{<:Any,<:Any,<:Array,<:Any,false}

⋮(x) = x
⋮(A::_LinIndexCPUArrayLike) = StridedView(A)
⋮(A::StructArray{T, N}) where {T,N} = StructArray{T, N}(⋮(StructArrays.components(A)))
⋮(tpl::Tuple) = map(⋮, tpl)
⋮(nt::NamedTuple) = map(⋮, nt)


@inline parallel_broadcast(f, args...) = broadcast(⋮(f), map(⋮, args)...)

@inline function parallel_broadcast!(f, result, args...)
    broadcast!(⋮(f), ⋮(result), map(⋮, args)...)
    return result
end


parallel_copyto!(A, B) = copyto!(⋮(A), ⋮(B))

parallel_copyto!(A::_CPUArrayLike, B::_NonLinIndexView) = _threads_parallel_copyto!(A, B)
parallel_copyto!(A::_NonLinIndexView, B::_CPUArrayLike) = _threads_parallel_copyto!(A, B)
parallel_copyto!(A::_NonLinIndexView, B::_NonLinIndexView) = _threads_parallel_copyto!(A, B)

function _threads_parallel_copyto!(A, B)
    idxs = eachindex(A)
    idxs == eachindex(B) || throw(ArgumentError("parallel_copyto! requires A and B to have exactly equal indices."))

    # ToDo: Tune heuristic for using multi-threading:
    if length(idxs) < 1000
        copyto!(A, B)
    else
        @inbounds Base.Threads.@threads for i in idxs
            A[i] = B[i]
        end
    end

    return A
end

function parallel_copyto!(A::StructArray{T,N}, B::StructArray{T,N}) where {T,N}
    map(parallel_copyto!, StructArrays.components(A), StructArrays.components(B))
    return A
end
