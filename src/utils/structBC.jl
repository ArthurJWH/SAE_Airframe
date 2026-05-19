abstract type AbstractBC end

struct NopBC <: AbstractBC end

struct FirstDerivativeBC{T} <: AbstractBC
    value::Float64
    index::T
end

struct SecondDerivativeBC{T} <: AbstractBC
    value::Float64
    index::T
end

struct ThirdDerivativeBC{T} <: AbstractBC
    value::Float64
    index::T
end