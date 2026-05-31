using GLMakie

using ..VLM: VLMMesh

"""
    plot_mesh(mesh::VLMMesh; min_extent = 5.0)

Interactive 3D visualization of a VLM mesh.

Arguments
- `mesh`: the VLM mesh to plot
- `min_extent`: minimum axis span in each direction
"""
function plot_mesh(mesh::VLMMesh; min_extent = 5.0)

    # Extract coordinates
    x = mesh.vertices[:, :, 1]
    y = mesh.vertices[:, :, 2]
    z = mesh.vertices[:, :, 3]

    # Helper to enforce a minimum axis range and include origin
    enforce_range(minv, maxv, minsize) = begin
        minv = min(minv, 0.0)
        maxv = max(maxv, 0.0)
        span = maxv - minv
        if span < minsize
            center = (minv + maxv) / 2
            half = minsize / 2
            return (center - half, center + half)
        end
        return (minv, maxv)
    end

    xlims = enforce_range(minimum(x), maximum(x), min_extent)
    ylims = enforce_range(minimum(y), maximum(y), min_extent)
    zlims = enforce_range(minimum(z), maximum(z), min_extent)

    # Create figure
    fig = GLMakie.Figure(
        size = (1000, 800)
    )

    # 3D axis
    ax = GLMakie.Axis3(
        fig[1,1],

        xlabel = "X",
        ylabel = "Y",
        zlabel = "Z",

        title = "VLM Mesh",

        # Equal scaling
        aspect = :data,

        # Nice perspective
        perspectiveness = 0.75
    )

    # Surface plot
    GLMakie.surface!(
        ax,
        x,
        y,
        z,

        shading = true,
        colormap = :viridis
    )

    if mesh.mirror_xz
        ylims = enforce_range(minimum(-y), maximum(y), min_extent)
        GLMakie.surface!(
            ax,
            x,
            -y,
            z,

            shading = true,
            colormap = :viridis
        )
        GLMakie.wireframe!(
        ax,
        x,
        -y,
        z,

        color = (:black, 0.4),
        linewidth = 1
    )
    end

    GLMakie.wireframe!(
        ax,
        x,
        y,
        z,

        color = (:black, 0.4),
        linewidth = 1
    )

    # Apply axis limits after creation
    GLMakie.xlims!(ax, xlims)
    GLMakie.ylims!(ax, ylims)
    GLMakie.zlims!(ax, zlims)

    # Display
    GLMakie.display(fig)

    return fig
end