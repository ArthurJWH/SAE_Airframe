module Setup

using Pkg

export setup

function setup(;
    dir = @__DIR__,
    reset = false,
    dependencies = ["Plots", "DataFrames"]
)

    dir = abspath(dir)

    if reset
        for file in ("Project.toml", "Manifest.toml")
            path = joinpath(dir, file)

            isfile(path) && rm(path)
        end
    end

    Pkg.activate(dir)

    Pkg.add(dependencies)

    Pkg.instantiate()

    println("Setup complete for: $dir")
end

end