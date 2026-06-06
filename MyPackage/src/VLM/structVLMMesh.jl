using ..Utils
using ..Geometry

struct VLMMesh
    vertices::Array{Float64,3}
    mirror_xz::Bool
end

function VLMMesh(plane::Plane, n_chord::Int, n_span::Int)
    meshes = [VLMMesh(surface, n_chord, n_span) for surface in plane.surfaces]
    return meshes
end

function VLMMesh(surface::Surface, n_chord::Int, n_span::Int)

    vertices = _generate_vertices(surface, n_chord, n_span, Val(surface.vertical))

    return VLMMesh(vertices, surface.mirror_xz)
end

@inline function _generate_vertices(surface::Surface, n_chord::Int, n_span::Int, ::Val{false})
    geom = _generate_geom(surface, n_chord, n_span)

    vertices = Array{Float64,3}(undef, n_chord+1, n_span+1, 3)

    @inbounds for i in eachindex(geom.x)
        xi = geom.x[i]
        spl_zi = geom.splines_z[i]
        for j in eachindex(geom.y)
            y = geom.y[j]
            chord_j = geom.chords[j]
            vertices[i, j, 1] = xi * chord_j
            vertices[i, j, 2] = y * geom.b
            vertices[i, j, 3] = spl_zi(y) * chord_j
        end
    end

    @inbounds for j in eachindex(geom.y)
        twist_j = geom.twists[j]
        sweep_length_j = geom.sweep_lengths[j]
        dihedral_length_j = geom.dihedral_lengths[j]

        _transform_section!(
            vertices,
            j,
            twist_j,
            geom.tw_center * geom.chords[j],
            sweep_length_j,
            geom.sw_center * (geom.root_chord - geom.chords[j]),
            dihedral_length_j,
            surface.pos
        )
    end

    return vertices
end

@inline function _generate_vertices(surface::Surface, n_chord::Int, n_span::Int, ::Val{true})
    geom = _generate_geom(surface, n_chord, n_span)

    vertices = Array{Float64,3}(undef, n_chord+1, n_span+1, 3)

    @inbounds for i in eachindex(geom.x)
        xi = geom.x[i]
        spl_zi = geom.splines_z[i]
        for j in eachindex(geom.y)
            y = geom.y[j]
            chord_j = geom.chords[j]
            vertices[i, j, 1] = xi * chord_j
            vertices[i, j, 3] = y * geom.b
            vertices[i, j, 2] = spl_zi(y) * chord_j
        end
    end

    @inbounds for j in eachindex(geom.y)
        twist_j = geom.twists[j]
        sweep_length_j = geom.sweep_lengths[j]
        dihedral_length_j = geom.dihedral_lengths[j]

        _transform_section_v!(
            vertices,
            j,
            twist_j,
            geom.tw_center * geom.chords[j],
            sweep_length_j,
            geom.sw_center * (geom.root_chord - geom.chords[j]),
            dihedral_length_j,
            surface.pos
        )
    end

    return vertices
end

@inline function _generate_geom(surface::Surface, n_chord::Int, n_span::Int)
    b = surface.b
    ys = surface.ys
    cambers = getfield.(surface.airfoils, :camber)
    chord = surface.chord
    twist = surface.twist
    tw_center = surface.tw_center
    sweep = surface.sweep
    sw_center = surface.sw_center
    dihedral = surface.dihedral
    
    x = (1 .- cos.(range(0, stop=pi, length=n_chord+1))) ./ 2
    y = (1 .- cos.(range(0, stop=pi, length=n_span+1))) ./ 2

    chords = chord.(y)
    
    # camber_grid = Matrix{Float64}(undef, length(x), length(ys))
    # @inbounds for j in eachindex(ys)
    #     camber_j = cambers[j]
    #     @simd for i in eachindex(x)
    #         camber_grid[i, j] = camber_j(x[i])
    #     end
    # end

    # splines_z = [
    #     LinearSpline(ys, @view camber_grid[i,:])
    #     for i in eachindex(x)
    # ]

    splines_z = Vector{LinearSpline}(undef, length(x))

    tmp = Vector{Float64}(undef, length(ys))

    @inbounds for i in eachindex(x)

        xi = x[i]

        for j in eachindex(ys)
            tmp[j] = cambers[j](xi)
        end

        splines_z[i] = LinearSpline(ys, copy(tmp))
    end

    sweep_integral = IntegrateGLQ(t -> tand(sweep(t)), n=5)
    dihedral_integral = IntegrateGLQ(t -> tand(dihedral(t)), n=5)
    sweep_lengths = Array{Float64}(undef, length(y))
    dihedral_lengths = Array{Float64}(undef, length(y))
    twists = Array{Float64}(undef, length(y))
    @inbounds for j in eachindex(y)
        yj = y[j]
        sweep_lengths[j] = b * sweep_integral(0.0, yj)
        dihedral_lengths[j] = b * dihedral_integral(0.0, yj)
        twists[j] = twist(yj)
    end

    return (
        b=b,
        ys=ys,
        tw_center=tw_center,
        sw_center=sw_center,
        x=x,
        y=y,
        chords=chords,
        root_chord=chords[1],
        splines_z=splines_z,
        sweep_lengths=sweep_lengths,
        dihedral_lengths=dihedral_lengths,
        twists=twists
    )
end

@inline function _transform_section!(
    vertices,
    j,
    twist,
    tw_center,
    sweep_length,
    offset,
    dihedral_length,
    pos
)

    c = cosd(twist)
    s = sind(twist)

    @inbounds @simd for i in axes(vertices,1)

        x0 = vertices[i,j,1]
        z0 = vertices[i,j,3]

        x = x0 - tw_center

        vertices[i,j,1] =
            x*c + z0*s + tw_center +
            offset + sweep_length + pos[1]

        vertices[i,j,2] += pos[2]

        vertices[i,j,3] =
            -x*s + z0*c + dihedral_length + pos[3]
    end

    return nothing
end

@inline function _transform_section_v!(
    vertices,
    j,
    twist,
    tw_center,
    sweep_length,
    offset,
    dihedral_length,
    pos
)

    c = cosd(twist)
    s = sind(twist)

    @inbounds @simd for i in axes(vertices,1)

        x0 = vertices[i,j,1]
        z0 = vertices[i,j,2]

        x = x0 - tw_center

        vertices[i,j,1] =
            x*c + z0*s + tw_center +
            offset + sweep_length + pos[1]

        vertices[i,j,2] = 
            -x*s + z0*c + dihedral_length + pos[2]

        vertices[i,j,3] += pos[3]
    end

    return nothing
end