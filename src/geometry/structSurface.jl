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
    MAC::Float64
    Swet::Float64

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