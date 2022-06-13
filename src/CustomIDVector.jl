

struct CustomIDArray{T, N} <: AbstractArray{T, N}
    data::Array{T, N}
    idx::Vector{Int}
end

const CustomIDVector{T} = CustomIDArray{T, 1}
CustomIDVector(data::Vector{T}, idx::Vector{Int}) where {T} = CustomIDVector{T}(data, idx)

const CustomIDMatrix{T} = CustomIDArray{T, 2}
CustomIDMatrix(data::Matrix{T}, idx::Vector{Int}) where {T} = CustomIDMatrix{T}(data, idx)

@inline size(cia::CustomIDArray) = size(cia.data)
@inline length(cia::CustomIDArray) = length(cia.data)


@inline function getindex(cia::CustomIDVector, i::Int)
    !(i in cia.idx) && throw(ArgumentError("No entry for contact with ID "*string(i)))
    idx::Int = Int(findfirst(ix -> ix == i, cia.idx))
    getindex(cia.data, idx)
end

@inline function getindex(cia::CustomIDArray{T,N}, I::Vararg{Int,N}) where {T, N}
    !all(in.(I, Ref(cia.idx))) && throw(ArgumentError("No entry for contact with ID "*string(I)))
    idx = broadcast(i -> findfirst(idx -> idx == i, cia.idx), I)
    getindex(cia.data, idx...)
end

@inline function push!(cia::CustomIDVector{T}, data::T, i::Int) where {T}
    push!(cia.data, data)
    push!(cia.idx, i)
end

@inline function setindex!(cia::CustomIDVector{T}, data::T, i::Int) where {T}
    if i in cia.idx
        idx::Int = Int(findfirst(idx -> idx == i, cia.idx))
        setindex!(cia.data, data, idx)
    else
        push!(cia, data, i)
    end
end

@inline function setindex!(cia::CustomIDArray{T, N}, data::T, I::Vararg{Int,N}) where {T, N}
    !all(in.(I, Ref(cia.idx))) && throw(ArgumentError("No entry for contact with ID "*string(I)))
    idx = broadcast(i -> findfirst(idx -> idx == i, cia.idx), I)
    setindex!(cia.data, data, idx...)
end

@inline function setindex!(cia::CustomIDArray{T, N}, data, I::Vararg{Int, N}) where {T, N}
    setindex!(cia, convert(T, data), I...)
end

@inline function Base.isequal(cia1::CIA1, cia2::CIA2) where {CIA1 <: CustomIDArray, CIA2 <: CustomIDArray}
    CIA1 == CIA2 && cia1.idx == cia2.idx && cia1.data == cia2.data
end

@inline function Base.:(==)(cia1::CIA1, cia2::CIA2) where {CIA1 <: CustomIDArray, CIA2 <: CustomIDArray}
    CIA1 == CIA2 && cia1.idx == cia2.idx && cia1.data == cia2.data
end

@inline iterate(cia::CustomIDArray, args...) = iterate(cia.data, args...)
@inline axes(cia::CustomIDVector) = (cia.idx,)
@inline size(cia::CustomIDArray) = size(cia.data)

@inline print(io::IO, cia::CustomIDArray) = print(io, cia.data)
@inline show(io::IO, cia::CustomIDArray) = print(io, cia)
@inline show(io::IO, ::MIME"text/plain", cia::CustomIDArray) = show(io, cia)
