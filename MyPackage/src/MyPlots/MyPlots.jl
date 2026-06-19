module MyPlots

using Plots

include("plot_airfoil.jl")
include("plot_mesh.jl")

export plot_airfoil
export plot_mesh

end
