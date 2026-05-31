using ..Coefficients

struct Plane
    surfaces::Vector{<:Surface}
    coeffs::Coeffs

    CG::Tuple{Float64, Float64, Float64}
    MTOW::Float64

    alpha_stall::Float64
    beta_stall::Float64

    SM::Float64
    alpha_trim::Float64
end

function Plane(surfaces::Vector{<:Surface})

    coeffs = Coeffs()
    CG = (0.0, 0.0, 0.0)
    MTOW = 0.0
    alpha_stall = 0.0
    beta_stall = 0.0
    SM = 0.0
    alpha_trim = 0.0

    return Plane(surfaces, coeffs, CG, MTOW, alpha_stall, beta_stall, SM, alpha_trim)
end