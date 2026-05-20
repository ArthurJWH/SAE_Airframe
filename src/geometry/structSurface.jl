include("../../assets/airfoils/Plain/Plain.jl")

using ..Utils

mutable struct Surface{chordF,twistF,sweepF,dihedralF}
    # y is the spanwise fraction coordinate
    # x is the chordwise fraction coordinate

    name::String
    mirror_xz::Bool
    pos::Tuple{Float64, Float64, Float64}
    rot::Tuple{Float64, Float64, Float64}

    b::Float64
    S::Float64
    AR::Float64
    MGC::Float64
    MAC::Float64
    # Swet::Float64 future drag build up

    ys::Vector{Float64}
    airfoils::Vector{String}

    # ys_mesh::Vector{Float64}

    chord::chordF
    twist::twistF
    tw_center::Float64
    sweep::sweepF
    sw_center::Float64
    dihedral::dihedralF
    
end

function Surface(;
    name::String="Surface",
    mirror_xz::Bool=true,
    pos::Tuple{Float64, Float64, Float64}=(0.0, 0.0, 0.0),
    rot::Tuple{Float64, Float64, Float64}=(0.0, 0.0, 0.0),
    b::Float64=1.0,
    ys::Vector{Float64}=[0.0, 1.0],
    airfoils::Vector{Airfoil}=[plain, plain],
    chord::chordF=y -> 1.0,
    twist::twistF=y -> 0.0,
    tw_center::Float64=0.25,
    sweep::sweepF=y -> 0.0,
    sw_center::Float64=0.25,
    dihedral::dihedralF=y -> 0.0
    ) where {chordF, twistF, sweepF, dihedralF}

    MGC = IntegrateGLQ(chord)(0, 1)
    S = MGC * b
    AR = b^2 / S
    MAC = 2 * IntegrateGLQ(x -> chord(x)^2)(0, 1) / MGC # Assume wing is symmetric about the centerline
    return Surface{chordF, twistF, sweepF, dihedralF}(name, mirror_xz, pos, rot,
                                                      b, S, AR, MGC, MAC, ys,
                                                      airfoils, chord, twist,
                                                      tw_center, sweep, sw_center,
                                                      dihedral)
end