using ..PlaneInfo

struct Plane
    surfaces::Vector{<:Surface}
    coeffs::Coeffs
    data::Data
end
