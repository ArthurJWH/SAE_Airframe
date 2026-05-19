module Utils

include("structBC.jl")
include("structSpline.jl")
include("structLSR.jl")
include("structInterpolate.jl")
include("structIntegrate.jl")

export AbstractBC, FirstDerivativeBC, SecondDerivativeBC, ThirdDerivativeBC, NaturalBC
export AbstractSpline, LinearSpline, QuadraticSpline, CubicSpline
export LSR
export Interpolate
export IntegrateGLQ

end