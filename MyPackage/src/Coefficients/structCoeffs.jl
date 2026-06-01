mutable struct Coeffs
    cl::Float64
    cd::Float64
    cm::Float64
end

Coeffs() = Coeffs(0.0, 0.0, 0.0)