# Internal functionality of SolidStateDetectors, not part of the public API:

const _CPUArrayLike = Union{Array, SubArray{<:Any,<:Any,<:Array}}
const _LinIndexCPUArrayLike = Union{Array, SubArray{<:Any,<:Any,<:Array,<:Any,true}}
const _NonLinIndexView = SubArray{<:Any,<:Any,<:Array,<:Any,false}


# Override for adapters that need special handling, to convert them to other adapters:
@inline pp_adapter(adapter::Adapter) where Adapter = adapter
@inline pp_adapter(::Type{Adapter}) where Adapter = Adapter
@inline pp_adapter(m::Module) = pp_module_adapter(Val(nameof(m)))

pp_module_adapter(Base.@nospecialize(module_name::Val)) = throw(ArgumentError("No default pp_adapter defined for module $(only(typeof(module_name).parameters))"))


const CurriedAdapt{Adapter} = Base.Fix1{typeof(Adapt.adapt), Adapter}

@inline adaptfunc(adapter::Adapter) where Adapter = Base.Fix1(Adapt.adapt, pp_adapter(adapter))
@inline adaptfunc(::Type{Adapter}) where Adapter = Base.Fix1(Adapt.adapt, pp_adapter(Adapter))

# ToDo (maybe): adaptfunct with multiple arguments, returning a Tuple?
# @inline adaptfunc(adapters::Vararg{Any,N}) where N = map(adaptfunc, adapters)


@inline adapted_call(f::F, af::CurriedAdapt, args...) where F = af(f)(map(af, args)...)
@inline adapted_call(f::F, adapter::A, args...) where {F,A} = adapted_call(f, adaptfunc(adapter), args...)
@inline adapted_call(f::F, ::Type{Adapter}, args...) where {F,Adapter} = adapted_call(f, adaptfunc(Adapter), args...)

@inline adapted_bcast(f::F, af::CurriedAdapt, args...) where F =  broadcast(af(f), map(af, args)...)
@inline adapted_bcast(f::F, adapter::A, args...) where {F,A} = adapted_bcast(f, adaptfunc(adapter), args...)
@inline adapted_bcast(f::F, ::Type{Adapter}, args...) where {F,Adapter} = adapted_bcast(f, adaptfunc(Adapter), args...)

@inline adapted_bcast!(f::F, af::CurriedAdapt, args...) where F = broadcast!(af(f), map(af, args)...)
@inline adapted_bcast!(f::F, adapter::A, args...) where {F,A} = adapted_bcast!(f, adaptfunc(adapter), args...)
@inline adapted_bcast!(f::F, ::Type{Adapter}, args...) where {F,Adapter} = adapted_bcast!(f, adaptfunc(Adapter), args...)


# For type adapters that need special handling:
struct PPTypeAdapter{T} end

const PPTypeAdaptFunc{T} = CurriedAdapt{PPTypeAdapter{T}}
Base.show(@nospecialize(io::IO), @nospecialize(f::PPTypeAdaptFunc{T})) where T = print(io, "adaptfunc($(nameof(T)))")
Base.show(@nospecialize(io::IO), ::MIME"text/plain", @nospecialize(f::PPTypeAdaptFunc{T})) where T = show(io, f)


@inline pp_adapter(::Type{Array}) = PPTypeAdapter{Array}()

Adapt.adapt_storage(::PPTypeAdapter{Array}, A::AbstractArray) = Array(A)



@inline pp_adapter(::Type{Strided.StridedView}) = PPTypeAdapter{Strided.StridedView}()
@inline pp_module_adapter(::Val{nameof(Strided)}) = pp_adapter(Strided.StridedView)
# Only adapt CPU arrays with linear indexing to StridedView:
Adapt.adapt_storage(::PPTypeAdapter{Strided.StridedView}, A::_LinIndexCPUArrayLike) = Strided.StridedView(A)


@inline pp_adapter(::Type{StrideArrays.StrideArray}) = PPTypeAdapter{StrideArrays.StrideArray}()
@inline pp_module_adapter(::Val{nameof(StrideArrays)}) = pp_adapter(StrideArrays.StrideArray)
# Only adapt CPU arrays with linear indexing to StrideArray:
Adapt.adapt_storage(::PPTypeAdapter{StrideArrays.StrideArray}, A::_LinIndexCPUArrayLike) = StrideArrays.StrideArray(A)


const ⋮ = adaptfunc(Strided.StridedView)
const ⋰ = adaptfunc(StrideArrays.StrideArray)


# ToDo: Add adapted_copyto! to replace parallel_copyto!?

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


# ToDo (maybe): Add: parallel_bcast and parallel_bcast like
# parallel_bcast(::Type{Strided.StridedView}, f, args...)
# parallel_bcast(::Type{StrideArrays.StrideArray}, f, args...)
# parallel_bcast(::Type{Task}, f, args...)
# parallel_bcast(::Type{Distributed.Worker}, f, args...)
# @inline is_parallelizing_arraytype(::Type{<:AbstractArray}) = Val(false)
# is_parallelizing_arraytype(T) = throw(ArgumentError("is_parallelizing_arraytype requires an array type as argument"))
# ...


# ToDo: Replace parallel_broadcast by adapted_bcast everywhere, then remove parallel_broadcast
@inline parallel_broadcast(f, args...) = adapted_bcast(f, ⋮, args...)
@inline parallel_broadcast!(f, args...) = adapted_bcast!(f, ⋮, args...)
