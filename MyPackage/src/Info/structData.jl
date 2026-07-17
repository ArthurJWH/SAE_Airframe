mutable struct Data
    CG::Tuple{Float64, Float64, Float64}
    SM::Float64
    MTOW::Float64
    alpha_stall::Float64
    beta_stall::Float64
    alpha_trim::Float64
end

function Data(;
    CG::Tuple{Float64, Float64, Float64}=(0.0, 0.0, 0.0),
    SM::Float64=0.0,
    MTOW::Float64=0.0,
    alpha_stall::Float64=0.0,
    beta_stall::Float64=0.0,
    alpha_trim::Float64=0.0,
)
    return Data(CG, SM, MTOW, alpha_stall, beta_stall, alpha_trim)
end
