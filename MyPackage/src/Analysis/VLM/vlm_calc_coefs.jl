function vlm_calc_coefs(
    L::Vector{Float64},
    D::Vector{Float64},
    M::Vector{Float64},
    Ml::Vector{Float64},
    N::Vector{Float64},
    rho::Float64,
    V_inf::Float64,
    S::Float64,
    MAC::Float64,
    b::Float64,
)
    Q = rho * V_inf^2 / 2

    CL = 0.0
    CD = 0.0

    CM = 0.0
    CMl = 0.0
    CN = 0.0

    for i in 1:length(L)
        CL += L[i] / (Q * S)
        CD += D[i] / (Q * S)

        CM += M[i] / (Q * S * MAC)
        CMl += Ml[i] / (Q * S * b)
        CN += N[i] / (Q * S * b)
    end

    return (CL, CD, CM, CMl, CN)
end
