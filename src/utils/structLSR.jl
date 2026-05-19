mutable struct LSR
    xs::Vector{<:AbstractFloat}
    fs::Vector{<:AbstractFloat}
    order::Int
    coeffs::Vector{<:AbstractFloat}
end

function LSR(coeffs::AbstractVector{<:AbstractFloat})
    return LSR([0.], [0.], length(coeffs)-1, coeffs)
end

function LSR(xs::AbstractVector{<:AbstractFloat}, fs::AbstractVector{<:AbstractFloat}, order::Int)
    @assert length(xs) == length(fs) "xs and fs must have same length"
    @assert issorted(xs) "xs must be sorted"
    lsr = LSR(xs, fs, order, zeros(order+1))
    polyfit!(lsr)
    return lsr
end

@inline function polyfit!(lsr::LSR)
    xs = lsr.xs
    M = zeros(length(xs), lsr.order+1)

    @views begin
        M[:, 1] .= 1
        for j in 2:lsr.order+1
            M[:, j] .= M[:, j-1] .* xs
        end

        @views lsr.coeffs .= M \ lsr.fs
    end
    return nothing
end

(lsr::LSR)(x::AbstractFloat) = evaluate(lsr, x)

@inline function evaluate(lsr::LSR, x::AbstractFloat)
    c = lsr.coeffs
    result = 0.0
    for i in length(c):-1:1
        result = result * x + c[i]
    end
    return result
end
