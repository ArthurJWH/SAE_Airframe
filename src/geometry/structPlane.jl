using ..Coefficients

struct Plane
    surfaces::Vector{Surface}
    coeffs::Coeffs

    CG::Tuple{Float64, Float64, Float64}
    MTOW::Float64

    alpha_stall::Float64
    beta_stall::Float64

    ME::Float64
    alpha_trim::Float64
end