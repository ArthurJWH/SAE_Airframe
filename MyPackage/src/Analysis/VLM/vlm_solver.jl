using LinearAlgebra
using StaticArrays
using Base.Threads

using ..Geometry

const Vec3 = SVector{3, Float64}

struct VLMSurface
    n_chord::Int
    n_span::Int
    mirror_xz::Bool
    range::UnitRange{Int}
end

struct VortexRing
    corners::NTuple{4, Vec3}
    colpt::Vec3
    normal::Vec3
    area::Float64
    surface_id::Int
end

struct GroundTransform
    h::Float64
    c2a::Float64
    s2a::Float64
    shift_x::Float64
    shift_z::Float64
    x_cg::Float64
    z_cg::Float64
    sa::Float64
end

function GroundTransform(rot::NTuple{2, Float64}, h::Float64, CG::Vec3)
    alpha = rot[1]
    c2a   = cosd(2 * alpha)
    s2a   = sind(2 * alpha)
    sa    = sind(alpha)
    return GroundTransform(
        h, c2a, s2a, 2 * h * s2a, 2 * h * c2a, CG[1], CG[3], sa
    )
end

struct VLMSetup
    initialized::Bool
    AIC::Matrix{Float64}
    n_panels::Int
    surfaces::Vector{VLMSurface}
    wake_map::Vector{Int}
    panel_rings::Vector{VortexRing}
    ground::Bool
    h::Float64
end

# function VLMSolver(
#     plane::Plane,
#     V_inf::Float64;
#     alpha::AbstractVector{<:Float64}=[0.0],
#     beta::AbstractVector{<:Float64}=[0.0],
#     ground::Bool=false,
#     h::Float64=0.0,
#     rho::Float64=1.225,
#     epsilon2::Float64=1e-10,
#     n_chord::Int=10,
#     n_span::Int=10,
#     wake_length::Float64=3.0
#     )
#     for b in beta
#         for a in alpha
#             VLMSolver(plane, V_inf, (a, b), ground=ground, h=h, rho=rho, epsilon2=epsilon2, n_chord=n_chord, n_span=n_span, wake_length=wake_length)
#         end
#     end
# end

function VLMSolver(
    plane::Plane,
    V_inf::Float64,
    rot::NTuple{2, Float64};
    ground::Bool=false,
    h::Float64=0.0,
    rho::Float64=1.225,
    epsilon2::Float64=1e-10,
    n_chordxspan::Vector{NTuple{2, Int}}=[(0, 0)],
)
    CG = Vec3(plane.data.CG...)

    meshes = VLMMesh(plane, n_chordxspan)

    rings, n_panels, n_surfaces, surfaces, wake_map = _gen_vortex_geom(meshes)

    if ground
        AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir, ring_img, dir_img = _assemble_sys(
            rings, n_panels, rot, wake_map, V_inf, h, CG, epsilon2
        )
    else
        AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir = _assemble_sys(
            rings, n_panels, rot, wake_map, V_inf, epsilon2
        )
    end

    setup = VLMSetup(
        true, AIC_rings, n_panels, surfaces, wake_map, rings, ground, h
    )

    gamma = AIC \ RHS

    if ground
        forces = _calc_forces(
            gamma,
            n_panels,
            n_surfaces,
            surfaces,
            rings,
            wake_map,
            rho,
            V_dir,
            V_inf,
            ring_img,
            dir_img,
            epsilon2,
        )
    else
        forces = _calc_forces(
            gamma,
            n_panels,
            n_surfaces,
            surfaces,
            rings,
            wake_map,
            rho,
            V_dir,
            V_inf,
            epsilon2,
        )
    end

    loads = _calc_loads(forces, rings, CG, n_surfaces, rot)

    # NEW: Trefftz-plane induced drag and total lift, appended as extra
    # return values.
    D_trefftz, L_trefftz = _calc_trefftz_loads(
        gamma, rings, surfaces, rho, V_dir, V_inf, rot, epsilon2
    )

    return (loads..., D_trefftz, L_trefftz)
end

function _gen_vortex_geom(meshes::AbstractVector{<:VLMMesh})
    n_meshes = length(meshes)
    surfaces = Vector{VLMSurface}(undef, n_meshes)
    rings    = VortexRing[]
    wake_map = Int[]

    sizehint!(rings, 1024) #TODO: improve sizehint
    sizehint!(wake_map, 128)

    i_start = 0
    for (surface_id, mesh) in enumerate(meshes)
        vertices = mesh.vertices
        mirror_xz = mesh.mirror_xz
        sz = size(vertices)
        # @assert sz[1] == 3 "First dimension must be 3 (x,y,z)"
        n_span = sz[2] - 1
        n_chord = sz[3] - 1

        _gen_mesh_rings!(
            rings, vertices, n_span, n_chord, surface_id, Val(mirror_xz)
        )

        i_end = i_start + (1 + mirror_xz) * n_span * n_chord

        surfaces[surface_id] = VLMSurface(
            n_chord, n_span, mirror_xz, (i_start + 1):i_end
        )
        append!(wake_map, (i_end - (1 + mirror_xz) * n_span + 1):i_end)

        i_start = i_end
    end

    n_panels = i_start
    return rings, n_panels, length(meshes), surfaces, wake_map
end

function _gen_mesh_rings!(
    rings::Vector{VortexRing},
    vertices::Array{Vec3, 3},
    n_span::Int,
    n_chord::Int,
    surface_id::Int,
    ::Val{false},
)::Nothing
    corners = Matrix{Vec3}(undef, n_span + 1, n_chord + 1)

    for i_chord in 1:n_chord
        @inbounds for i_span in 1:(n_span + 1)
            le = Vec3(
                vertices[1, i_span, i_chord],
                vertices[2, i_span, i_chord],
                vertices[3, i_span, i_chord],
            )
            te = Vec3(
                vertices[1, i_span, i_chord + 1],
                vertices[2, i_span, i_chord + 1],
                vertices[3, i_span, i_chord + 1],
            )
            corners[i_span, i_chord] = le + 0.25 * (te - le)
        end
    end
    @inbounds for i_span in 1:(n_span + 1)
        corners[i_span, n_chord + 1] = Vec3(
            vertices[1, i_span, n_chord + 1],
            vertices[2, i_span, n_chord + 1],
            vertices[3, i_span, n_chord + 1],
        )
    end

    @inbounds for i_chord in 1:n_chord
        for i_span in 1:n_span
            A = corners[i_span, i_chord]
            B = corners[i_span + 1, i_chord]
            C = corners[i_span + 1, i_chord + 1]
            D = corners[i_span, i_chord + 1]

            colpt = 0.25 * (A + B + C + D)
            n, area = _panel_normal_area(A, B, C, D)

            rings[(i_chord - 1) * n_span + i_span] = VortexRing(
                (A, B, C, D), colpt, n, area, surface_id
            )
        end
    end

    return nothing
end

function _gen_mesh_rings!(
    rings::Vector{VortexRing},
    vertices::Array{Vec3, 3},
    n_span::Int,
    n_chord::Int,
    surface_id::Int,
    ::Val{true},
)::Nothing
    corners = Matrix{Vec3}(undef, n_span + 1, n_chord + 1)
    rings_helper = Vector{VortexRing}(undef, 2 * n_span * n_chord)

    for i_chord in 1:n_chord
        @inbounds for i_span in 1:(n_span + 1)
            le = Vec3(
                vertices[1, i_span, i_chord],
                vertices[2, i_span, i_chord],
                vertices[3, i_span, i_chord],
            )
            te = Vec3(
                vertices[1, i_span, i_chord + 1],
                vertices[2, i_span, i_chord + 1],
                vertices[3, i_span, i_chord + 1],
            )
            corners[i_span, i_chord] = le + 0.25 * (te - le)
        end
    end
    @inbounds for i_span in 1:(n_span + 1)
        corners[i_span, n_chord + 1] = Vec3(
            vertices[1, i_span, n_chord + 1],
            vertices[2, i_span, n_chord + 1],
            vertices[3, i_span, n_chord + 1],
        )
    end

    @inbounds for i_chord in 1:n_chord
        for i_span in 1:n_span
            A = corners[i_span, i_chord]
            B = corners[i_span + 1, i_chord]
            C = corners[i_span + 1, i_chord + 1]
            D = corners[i_span, i_chord + 1]
            colpt = 0.25 * (A + B + C + D)
            n, area = _panel_normal_area(A, B, C, D)

            rings_helper[(i_chord - 1) * 2 * n_span + i_span] = VortexRing(
                (A, B, C, D), colpt, n, area, surface_id
            )

            A_m = (A[1], -A[2], A[3])
            B_m = (B[1], -B[2], B[3])
            C_m = (C[1], -C[2], C[3])
            D_m = (D[1], -D[2], D[3])
            colpt_m = Vec3(colpt[1], -colpt[2], colpt[3])
            n_m = Vec3(n[1], -n[2], n[3])

            rings_helper[(i_chord - 1) * 2 * n_span + i_span + n_span] = VortexRing(
                (A_m, B_m, C_m, D_m), colpt_m, n_m, area, surface_id
            )
        end
    end

    append!(rings, rings_helper)

    return nothing
end

@inline function _panel_normal_area(
    A::Vec3, B::Vec3, C::Vec3, D::Vec3
)::Tuple{Vec3, Float64}
    d1 = C - A   # diagonal 1
    d2 = D - B   # diagonal 2
    cr = cross(d1, d2)
    normal = norm(cr)
    n = cr / normal
    area = 0.5 * normal
    return n, area
end

function _assemble_sys(
    rings::Vector{VortexRing},
    n_panels::Int,
    rot::NTuple{2, Float64},
    wake_map::Vector{Int},
    V_inf::Float64,
    epsilon2::Float64,
)
    V_dir = _wind_dir(rot)
    inv_wake_map = _build_inv_wake_map(wake_map, n_panels)

    AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz = _assemble_AIC(
        rings, n_panels, inv_wake_map, V_dir, epsilon2
    )

    RHS = _assemble_RHS(rings, n_panels, V_dir, V_inf)

    return AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir
end

function _assemble_sys(
    rings::Vector{VortexRing},
    n_panels::Int,
    rot::NTuple{2, Float64},
    wake_map::Vector{Int},
    V_inf::Float64,
    h::Float64,
    CG::Vec3,
    epsilon2::Float64,
)
    V_dir = _wind_dir(rot)
    inv_wake_map = _build_inv_wake_map(wake_map, n_panels)

    gt       = GroundTransform(rot, h, CG)
    ring_img = _precompute_image_corners(rings, gt)
    dir_img  = _apply_ground_transform_dir(V_dir, gt)

    AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz = _assemble_AIC(
        rings, n_panels, inv_wake_map, V_dir, ring_img, dir_img, epsilon2
    )

    RHS = _assemble_RHS(rings, n_panels, V_dir, V_inf)

    return AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir, ring_img, dir_img
end

@inline function _wind_dir(rot::NTuple{2, Float64})::Vec3
    alpha, beta = rot
    ca, sa = cosd(alpha), sind(alpha)
    cb, sb = cosd(beta), sind(beta)
    return Vec3(ca * cb, -sb, sa * cb)
end

function _build_inv_wake_map(wake_map::Vector{Int}, n_panels::Int)::Vector{Int}
    inv = zeros(Int, n_panels)
    @inbounds for k in eachindex(wake_map)
        inv[wake_map[k]] = k
    end
    return inv
end

function _precompute_image_corners(
    panels::Vector{VortexRing}, gt::GroundTransform
)::Vector{NTuple{4, Vec3}}
    n   = length(panels)
    img = Vector{NTuple{4, Vec3}}(undef, n)
    @inbounds for j in 1:n
        img[j] = _apply_ground_transform(panels[j].corners, gt)
    end
    return img
end

@inline function _apply_ground_transform(
    corners::NTuple{4, Vec3}, gt::GroundTransform
)::NTuple{4, Vec3}
    ntuple(Val(4)) do i
        x, y, z = corners[i]
        x_ref, z_ref = x - gt.x_cg, z - gt.z_cg
        x_t = gt.c2a * x_ref - gt.s2a * z_ref + gt.shift_x + gt.x_cg
        z_t = gt.s2a * x_ref + gt.c2a * z_ref - gt.shift_z + gt.z_cg
        # @assert z_t < z "Mirrored point is above the original point, check h and CG values."
        @assert gt.h > gt.sa * x_ref "Plane intersects the ground, increase h or adjust CG."
        Vec3(x_t, y, z_t)
    end
end

@inline function _apply_ground_transform_dir(
    dir::Vec3, gt::GroundTransform
)::Vec3
    x, y, z = dir
    x_t = gt.c2a * x - gt.s2a * z
    z_t = gt.s2a * x + gt.c2a * z
    return Vec3(x_t, y, z_t)
end

function _assemble_AIC(
    rings::Vector{VortexRing},
    n_panels::Int,
    inv_wake_map::Vector{Int},
    V_dir::Vec3,
    epsilon2::Float64,
)
    AIC = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_ring = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vx = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vy = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vz = Matrix{Float64}(undef, n_panels, n_panels)

    # Parallelise over rows (collocation points)
    @threads for i in 1:n_panels
        ri = rings[i]
        colpt = ri.colpt
        ni = ri.normal

        @inbounds for j in 1:n_panels
            corners = rings[j].corners
            vel = _ring_induced_v(colpt, corners, 1.0, epsilon2)

            AIC_vx[i, j] = vel[1]
            AIC_vy[i, j] = vel[2]
            AIC_vz[i, j] = vel[3]

            k = inv_wake_map[j]
            if k > 0
                vel += _horseshoe_induced_v(
                    colpt, corners[4], corners[3], V_dir, 1.0, epsilon2
                )
            end

            AIC[i, j] = dot(ni, vel)
        end
    end

    return AIC, AIC_ring, AIC_vx, AIC_vy, AIC_vz
end

function _assemble_AIC(
    rings::Vector{VortexRing},
    n_panels::Int,
    inv_wake_map::Vector{Int},
    V_dir::Vec3,
    ring_img::Vector{NTuple{4, Vec3}},
    dir_img::Vec3,
    epsilon2::Float64,
)
    AIC = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_ring = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vx = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vy = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vz = Matrix{Float64}(undef, n_panels, n_panels)

    # Parallelise over rows (collocation points)
    @threads for i in 1:n_panels
        ri = rings[i]
        colpt = ri.colpt
        ni = ri.normal

        @inbounds for j in 1:n_panels
            rj = rings[j]
            corners = rj.corners
            img_corners = ring_img[j]

            vel =
                _ring_induced_v(colpt, corners, 1.0, epsilon2) -
                _ring_induced_v(colpt, img_corners, 1.0, epsilon2)

            AIC_vx[i, j] = vel[1]
            AIC_vy[i, j] = vel[2]
            AIC_vz[i, j] = vel[3]

            AIC_ring[i, j] = dot(ni, vel)

            k = inv_wake_map[j]
            if k > 0
                vel +=
                    _horseshoe_induced_v(
                        colpt, corners[4], corners[3], V_dir, 1.0, epsilon2
                    ) - _horseshoe_induced_v(
                        colpt,
                        img_corners[4],
                        img_corners[3],
                        dir_img,
                        1.0,
                        epsilon2,
                    )
            end

            AIC[i, j] = dot(ni, vel)
        end
    end

    return AIC, AIC_ring, AIC_vx, AIC_vy, AIC_vz
end

@inline function _segment_induced_v(
    P::Vec3, A::Vec3, B::Vec3, Gamma::Float64, epsilon2::Float64
)::Vec3
    r1 = A - P
    r2 = B - P
    r0 = B - A
    r0_2 = dot(r0, r0)

    # degenerate (zero-length) segment guard, independent of P.
    if r0_2 < 1e-15 # or 1e-28
        return Vec3(0.0, 0.0, 0.0)
    end

    # epsilon2 is a dimensionless core-radius fraction of the
    # local segment length^2 ([core_len2] = m^2), so the regularization
    # scales with panel size instead of being a fixed, geometry-independent
    # absolute number.
    core2 = epsilon2 * r0_2

    cr = cross(r1, r2)
    cr2 = dot(cr, cr) + core2 * core2

    # norm1 = norm(r1)
    # norm2 = norm(r2)

    # if cr2 == 0.0 || norm1 < 1e-14 || norm2 < 1e-14
    #     return Vec3(0.0, 0.0, 0.0)
    # end

    norm1 = sqrt(dot(r1, r1) + core2)
    norm2 = sqrt(dot(r2, r2) + core2)

    factor =
        Gamma / (4 * pi * cr2) * (dot(r0, r1) / norm1 - dot(r0, r2) / norm2)
    return factor * cr
end

@inline function _semi_infinite_induced_v(
    P::Vec3,
    A::Vec3,
    dir::Vec3,
    Gamma::Float64,
    ref_len2::Float64,
    epsilon2::Float64,
)::Vec3
    if ref_len2 < 1e-15 # or 1e-28
        return Vec3(0.0, 0.0, 0.0)
    end

    r1 = A - P
    r1_2 = dot(r1, r1)
    core2 = epsilon2 * ref_len2

    cr = cross(r1, dir)
    cr2 = dot(cr, cr) + core2 * core2

    norm1 = sqrt(r1_2 + core2)
    # if norm1 < 1e-14
    #     return Vec3(0.0, 0.0, 0.0)
    # end
    cos_th = dot(r1, dir) / norm1

    factor = Gamma / (4 * pi * cr2) * (cos_th - 1.0)
    return factor * cr
end

function _ring_induced_v(
    P::Vec3,
    A::Vec3,
    B::Vec3,
    C::Vec3,
    D::Vec3,
    Gamma::Float64,
    epsilon2::Float64,
)::Vec3
    v1 = _segment_induced_v(P, A, B, Gamma, epsilon2)
    v2 = _segment_induced_v(P, B, C, Gamma, epsilon2)
    v3 = _segment_induced_v(P, C, D, Gamma, epsilon2)
    v4 = _segment_induced_v(P, D, A, Gamma, epsilon2)

    return v1 + v2 + v3 + v4
end

@inline function _ring_induced_v(
    P::Vec3, corners::NTuple{4, Vec3}, Gamma::Float64, epsilon2::Float64
)::Vec3
    return _ring_induced_v(
        P, corners[1], corners[2], corners[3], corners[4], Gamma, epsilon2
    )
end

@inline function _horseshoe_induced_v(
    P::Vec3, A::Vec3, B::Vec3, dir::Vec3, Gamma::Float64, epsilon2::Float64
)::Vec3
    ref_len2 = dot(B - A, B - A)
    return _segment_induced_v(P, A, B, Gamma, epsilon2) +
           _semi_infinite_induced_v(P, B, dir, Gamma, ref_len2, epsilon2) -
           _semi_infinite_induced_v(P, A, dir, Gamma, ref_len2, epsilon2)
end

function _assemble_RHS(
    rings::Vector{VortexRing}, n_panels::Int, V_dir::Vec3, V_inf::Float64
)::Vector{Float64}
    RHS = Vector{Float64}(undef, n_panels)

    @threads for i in 1:n_panels
        RHS[i] = -dot(rings[i].normal, V_dir * V_inf)
    end

    return RHS
end

function _calc_forces(
    gamma::Vector{Float64},
    n_panels::Int,
    n_surfaces::Int,
    surfaces::Vector{VLMSurface},
    rings::Vector{VortexRing},
    wake_map::Vector{Int},
    rho::Float64,
    V_dir::Vec3,
    V_inf::Float64,
    epsilon2::Float64,
)
    n_wakes = length(wake_map)
    V_free = V_dir * V_inf
    forces = Vector{Matrix{Vec3}}(undef, n_surfaces)

    # Compute induced velocity at every collocation point
    # (sum of contributions from all lifting panels and all wake rings)
    V_ind = Vector{Vec3}(undef, n_panels)

    @threads for i in 1:n_panels
        ri = rings[i]
        colpt = ri.colpt
        vel = Vec3(0.0, 0.0, 0.0)

        @inbounds for j in 1:n_panels
            vel += _ring_induced_v(colpt, rings[j].corners, gamma[j], epsilon2)
        end

        @inbounds for j in 1:n_wakes
            p = wake_map[j]
            pc = rings[p].corners
            vel += _horseshoe_induced_v(
                colpt, pc[4], pc[3], V_dir, gamma[p], epsilon2
            )
        end

        V_ind[i] = vel
    end

    for i in 1:n_surfaces
        surface = surfaces[i]
        n_span = surface.n_span
        n_chord = surface.n_chord
        mirror_xz = surface.mirror_xz
        range = surface.range
        start = range.start

        forces[i] = Matrix{Vec3}(undef, (mirror_xz + 1) * n_span, n_chord)

        for j in range
            rj = rings[j]
            span_vec = rj.corners[2] - rj.corners[1]
            chord_vec = rj.corners[4] - rj.corners[1]

            delta_gamma_c, delta_gamma_s = _calc_delta_gamma(
                gamma, j, start, n_span, mirror_xz
            )

            V_total = V_free + V_ind[j]
            vector = delta_gamma_c * span_vec + delta_gamma_s * chord_vec

            i_span, i_chord = _j1dto2d(j, start, n_span, mirror_xz)
            forces[i][i_span, i_chord] = rho * cross(V_total, vector)
        end
    end

    return forces
end

function _calc_forces(
    gamma::Vector{Float64},
    n_panels::Int,
    n_surfaces::Int,
    surfaces::Vector{VLMSurface},
    rings::Vector{VortexRing},
    wake_map::Vector{Int},
    rho::Float64,
    V_dir::Vec3,
    V_inf::Float64,
    ring_img::Vector{NTuple{4, Vec3}},
    dir_img::Vec3,
    epsilon2::Float64,
)
    n_wakes = length(wake_map)
    V_free = V_dir * V_inf
    forces = Vector{Matrix{Vec3}}(undef, n_surfaces)

    # Compute induced velocity at every collocation point
    # (sum of contributions from all lifting panels and all wake rings)
    V_ind = Vector{Vec3}(undef, n_panels)

    @threads for i in 1:n_panels
        ri = rings[i]
        colpt = ri.colpt
        vel = Vec3(0.0, 0.0, 0.0)

        @inbounds for j in 1:n_panels
            vel +=
                _ring_induced_v(colpt, rings[j].corners, gamma[j], epsilon2) -
                _ring_induced_v(colpt, ring_img[j], gamma[j], epsilon2)
        end

        # CHANGED: real horseshoe from the TE panel's own corners + V_dir;
        # image horseshoe from ring_img's mirrored TE corners + dir_img.
        @inbounds for j in 1:n_wakes
            p = wake_map[j]
            pc = rings[p].corners
            pic = ring_img[p]
            vel +=
                _horseshoe_induced_v(
                    colpt, pc[4], pc[3], V_dir, gamma[p], epsilon2
                ) - _horseshoe_induced_v(
                    colpt, pic[4], pic[3], dir_img, gamma[p], epsilon2
                )
        end

        V_ind[i] = vel
    end

    for i in 1:n_surfaces
        surface = surfaces[i]
        n_span = surface.n_span
        n_chord = surface.n_chord
        mirror_xz = surface.mirror_xz
        range = surface.range
        start = range.start

        forces[i] = Matrix{Vec3}(undef, (mirror_xz + 1) * n_span, n_chord)

        for j in range
            rj = rings[j]
            span_vec = rj.corners[2] - rj.corners[1]
            chord_vec = rj.corners[4] - rj.corners[1]

            delta_gamma_c, delta_gamma_s = _calc_delta_gamma(
                gamma, j, start, n_span, Val(mirror_xz)
            )

            V_total = V_free + V_ind[j]
            vector = delta_gamma_c * span_vec + delta_gamma_s * chord_vec

            i_span, i_chord = _j1dto2d(j, start, n_span, Val(mirror_xz))
            forces[i][i_span, i_chord] = rho * cross(V_total, vector)
        end
    end

    return forces
end

function _calc_delta_gamma(
    gamma::Vector{Float64}, j::Int, start::Int, n_span::Int, ::Val{false}
)
    delta_gamma_s = gamma[j] - ((j - start) % n_span == 0 ? 0.0 : gamma[j - 1])
    delta_gamma_c = gamma[j] - ((j - start) < n_span ? 0.0 : gamma[j - n_span])
    return delta_gamma_c, delta_gamma_s
end

function _calc_delta_gamma(
    gamma::Vector{Float64}, j::Int, start::Int, n_span::Int, ::Val{true}
)
    delta_gamma_s =
        gamma[j] - ((j - start) % (2 * n_span) == 0 ? 0.0 : gamma[j - 1])
    delta_gamma_c =
        gamma[j] - ((j - start) < (2 * n_span) ? 0.0 : gamma[j - 2 * n_span])
    return delta_gamma_c, delta_gamma_s
end

function _j1dto2d(j::Int, start::Int, n_span::Int, ::Val{false})
    i_span = (j - start) % n_span + 1
    i_chord = (j - start)//n_span + 1

    return i_span, i_chord
end

function _j1dto2d(j::Int, start::Int, n_span::Int, ::Val{true})
    i_span = (j - start) % (2 * n_span) + 1
    i_chord = (j - start)//(2 * n_span) + 1

    return i_span, i_chord
end

function _calc_loads(
    forces::Vector{Matrix{Vec3}},
    rings::Vector{VortexRing},
    CG::Vec3,
    n_surfaces::Int,
    rot::NTuple{2, Float64},
)
    alpha = rot[1]
    beta = rot[2]
    sa = sind(alpha)
    ca = cosd(alpha)
    sb = sind(beta)
    cb = cosd(beta)

    FX = fill(0.0, n_surfaces) # Force in x direction
    FY = fill(0.0, n_surfaces) # Force in y direction
    FZ = fill(0.0, n_surfaces) # Force in z direction

    FX_dist = Vector{Vector{Float64}}(undef, n_surfaces)
    FY_dist = Vector{Vector{Float64}}(undef, n_surfaces)
    FZ_dist = Vector{Vector{Float64}}(undef, n_surfaces)

    L = fill(0.0, n_surfaces) # Lift
    D = fill(0.0, n_surfaces) # Drag

    L_dist = Vector{Vector{Float64}}(undef, n_surfaces)
    D_dist = Vector{Vector{Float64}}(undef, n_surfaces)

    M = fill(0.0, n_surfaces) # Pitching moment
    Ml = fill(0.0, n_surfaces) # Roll moment
    N = fill(0.0, n_surfaces) # Yaw moment

    M_dist = Vector{Vector{Float64}}(undef, n_surfaces)
    Ml_dist = Vector{Vector{Float64}}(undef, n_surfaces)
    N_dist = Vector{Vector{Float64}}(undef, n_surfaces)

    i1d = 1

    for i in 1:n_surfaces
        surface_forces = forces[i]
        n_span, n_chord = size(surface_forces)

        FX_dist[i] = fill(0.0, n_span)
        FY_dist[i] = fill(0.0, n_span)
        FZ_dist[i] = fill(0.0, n_span)

        L_dist[i] = fill(0.0, n_span)
        D_dist[i] = fill(0.0, n_span)

        M_dist[i] = fill(0.0, n_span)
        Ml_dist[i] = fill(0.0, n_span)
        N_dist[i] = fill(0.0, n_span)

        for j in 1:n_span
            for k in 1:n_chord
                fx = surface_forces[j, k][1]
                fy = surface_forces[j, k][2]
                fz = surface_forces[j, k][3]

                FX[i] += fx
                FY[i] += fy
                FZ[i] += fz

                FX_dist[i][j] += fx
                FY_dist[i][j] += fy
                FZ_dist[i][j] += fz

                lift = -fx * sa + fz * ca
                drag = fx * ca * cb - fy * sb + fz * sa * cb

                L[i] += lift
                D[i] += drag

                L_dist[i][j] += lift
                D_dist[i][j] += drag

                rj = rings[i1d]
                i1d += 1
                moment_arm = rj.colpt - CG

                (ml, m, n) = cross(moment_arm, surface_forces[j, k])

                M[i] += m
                Ml[i] += ml
                N[i] += n

                M_dist[i][j] += m
                Ml_dist[i][j] += ml
                N_dist[i][j] += n
            end
        end
    end

    return (
        FX,
        FX_dist,
        FY,
        FY_dist,
        FZ,
        FZ_dist,
        L,
        L_dist,
        D,
        D_dist,
        M,
        M_dist,
        Ml,
        Ml_dist,
        N,
        N_dist,
    )
end

# NEW: induced drag AND total lift via Trefftz-plane (far-field) analysis,
# as a cross-check against the near-field Kutta-Joukowski summation in
# _calc_forces. Far downstream, all trailing legs are parallel to V_dir, so
# the problem reduces to a 2D vortex sheet in the plane perpendicular to
# V_dir; this avoids the near-field method's sensitivity to local induced-
# velocity noise close to the panels/tip vortices.
#
# Drag: (1) build the trailing-vortex sheet as a set of semi-infinite
# filaments at each spanwise TE panel boundary, with strength equal to the
# local circulation jump (root and tip edges get the panel's full local
# gamma, since nothing cancels them there); (2) at each TE panel's own
# spanwise station, evaluate the far-field (Trefftz-plane) cross-flow
# velocity induced by that whole filament sheet, using the exact far-field
# limit of _semi_infinite_induced_v; (3) apply the same Kutta-Joukowski
# cross-product convention already used in _calc_forces, but with the
# freestream term dropped (only the induced velocity does work along
# V_dir, i.e. produces drag).
#
# Lift: the classical Trefftz-plane/Kutta-Joukowski result is that
# total lift depends ONLY on the total bound circulation and V_free, not
# on the local downwash w computed for drag above -- so it's obtained with
# the same cross-product form but using V_free in place of w, then
# extracted with the exact same lift-axis projection _calc_loads uses
# (lift = -Fx*sin(alpha) + Fz*cos(alpha)), so L_trefftz is directly
# comparable to the near-field lift returned by _calc_loads.
#
# NOTE: this assumes the global (y,z) cross-flow plane is a reasonable
# proxy for the true Trefftz plane, i.e. small-to-moderate alpha/beta --
# consistent with the flat-wake, linear assumptions.

function _calc_trefftz_loads(
    gamma::Vector{Float64},
    rings::Vector{VortexRing},
    surfaces::Vector{VLMSurface},
    rho::Float64,
    V_dir::Vec3,
    V_inf::Float64,
    rot::NTuple{2, Float64},
    epsilon2::Float64,
)::Tuple{Float64, Float64}
    fil_pos = Vec3[]
    fil_strength = Float64[]
    fil_ref2 = Float64[]

    for surf in surfaces
        n_span = surf.n_span
        te_start = surf.range.stop - (1 + surf.mirror_xz) * n_span + 1
        for half in 0:(surf.mirror_xz)
            base = te_start + half * n_span
            prev_g = 0.0
            for k in 1:n_span
                j = base + k - 1
                g = gamma[j]
                rj = rings[j]
                edge2 = dot(
                    rj.corners[3] - rj.corners[4], rj.corners[3] - rj.corners[4]
                )
                push!(fil_pos, rj.corners[4])  # root-ward edge
                push!(fil_strength, g - prev_g)
                push!(fil_ref2, edge2)
                prev_g = g
                if k == n_span
                    push!(fil_pos, rj.corners[3])  # tip-closing edge
                    push!(fil_strength, -g)
                    push!(fil_ref2, edge2)
                end
            end
        end
    end
    n_fil = length(fil_pos)

    V_free = V_dir * V_inf
    sa, ca = sind(rot[1]), cosd(rot[1])

    D = 0.0
    F_trefftz = Vec3(0.0, 0.0, 0.0)
    for surf in surfaces
        n_span = surf.n_span
        te_start = surf.range.stop - (1 + surf.mirror_xz) * n_span + 1
        for half in 0:(surf.mirror_xz)
            base = te_start + half * n_span
            for k in 1:n_span
                j = base + k - 1
                rj = rings[j]
                span_vec = rj.corners[3] - rj.corners[4]
                midpt = 0.5 * (rj.corners[3] + rj.corners[4])

                w = Vec3(0.0, 0.0, 0.0)
                @inbounds for m in 1:n_fil
                    r1 = fil_pos[m] - midpt
                    r1 -= dot(r1, V_dir) * V_dir
                    cr = cross(r1, V_dir)
                    h2 = dot(cr, cr) + epsilon2 * fil_ref2[m]
                    w += (-fil_strength[m] / (4 * pi * h2)) * cr
                end

                D += dot(rho * cross(w, gamma[j] * span_vec), V_dir)
                F_trefftz += rho * cross(V_free, gamma[j] * span_vec)
            end
        end
    end

    L = -F_trefftz[1] * sa + F_trefftz[3] * ca

    return D, L
end
