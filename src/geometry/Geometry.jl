module Geometry

include("structAirfoil.jl")
include("structSurface.jl")
include("structPlane.jl")

export Airfoil, calc_surfaces, calc_camber
export Surface
export Plane

end