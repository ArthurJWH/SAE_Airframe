include("src/Setup.jl")

using .Setup

setup(dir = @__DIR__, reset = false, dependencies = ["Plots"])