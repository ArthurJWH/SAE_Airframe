using ..Utils
using ..Geometry

struct VLMMesh
    vertices::Array{Float64,3}
end

function VLMMesh(surface::Surface, n_chord::Int, n_span::Int)

    b = surface.b
    ys = surface.ys
    airfoils = surface.airfoils
    cambers = tuple([airfoil.camber for airfoil in airfoils]...)
    chord = surface.chord
    twist = surface.twist
    tw_center = surface.tw_center
    sweep = surface.sweep
    sw_center = surface.sw_center
    dihedral = surface.dihedral

    vertices = Array{Float64,3}(undef, n_chord+1, n_span+1, 3)
    x = (1 .- cos.(range(0, stop=pi, length=n_chord+1))) ./ 2
    y = (1 .- cos.(range(0, stop=pi, length=n_span+1))) ./ 2
    # z = interpolate the camber line for each airfoil at the given x coordinates
    # for each x, calculate z for each section and spline and calculate z for each section with y
    chords = chord.(y)

    camber_grid = Matrix{Float64}(undef, length(x), length(ys))

    @inbounds for j in eachindex(ys)
        camber_j = cambers[j]
        for i in eachindex(x)
            camber_grid[i, j] = camber_j(x[i])
        end
    end

    @inbounds for i in eachindex(x)
        spl_z = LinearSpline(ys, view(camber_grid, i, :))
        for j in eachindex(y)
            chord_j = chords[j]
            vertices[i,j,1] = x[i] * chord_j
            vertices[i,j,2] = y[j] * b
            vertices[i,j,3] = spl_z(y[j]) * chord_j
        end
    end

    sweep_length(y) = IntegrateGLQ(sweep, n=5)(0., y)
    dihedral_length(y) = IntegrateGLQ(dihedral, n=5)(0., y)

    @inbounds for j in eachindex(y)

        # apply twist, sweep, and dihedral to each vertex in the section
        twist_j = twist(y[j])
        sweep_length_j = sweep_length(y[j])
        dihedral_length_j = dihedral_length(y[j])

        section_vertices = @views vertices[:, j, :]
        # apply twist
        rotate_section!(section_vertices, twist_j, tw_center * chords[j])
        # apply sweep
        sweep_section!(section_vertices, sweep_length_j, sw_center * (chords[1] -chords[j]))
        # apply dihedral
        dihedral_section!(section_vertices, dihedral_length_j)

    end

    return VLMMesh(vertices)
end

@inline function rotate_section!(section_vertices, twist, tw_center)

    c = cosd(twist)
    s = sind(twist)

    @inbounds for i in axes(section_vertices,1)

        x = section_vertices[i,1] - tw_center
        z = section_vertices[i,3]

        section_vertices[i,1] = x*c + z*s + tw_center
        section_vertices[i,3] = -x*s + z*c

    end

    return nothing
end

@inline function sweep_section!(section_vertices, sweep_length, offset)
    section_vertices[:, 1] .+= offset + sweep_length
    return nothing
end

@inline function dihedral_section!(section_vertices, dihedral_length)
    section_vertices[:, 3] .+= dihedral_length
    return nothing
end