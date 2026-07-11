mutable struct Interpolate
    xs::Vector{Float64}
    fs::Vector{Float64}
    coeffs::Vector{Float64}
end

function Interpolate(coeffs::AbstractVector{<:AbstractFloat})
    return Interpolate([0.0], [0.0], coeffs)
end

function Interpolate(
    xs::AbstractVector{<:AbstractFloat}, fs::AbstractVector{<:AbstractFloat}
)
    @assert length(xs) == length(fs) "xs and fs must have same length"
    @assert issorted(xs) "xs must be sorted"
    n = length(xs)
    return LSR(xs, fs, n - 1)
end

(interp::Interpolate)(x::AbstractFloat) = evaluate(interp, x)

@inline function evaluate(interp::Interpolate, x::AbstractFloat)
    c = interp.coeffs
    result = 0.0
    for i in length(c):-1:1
        result = result * x + c[i]
    end
    return result
end
