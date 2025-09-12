ENV["GKSwstype"] = "100"

using Pkg

Pkg.add("DrWatson")

Pkg.activate(".")

if Pkg.TOML.parsefile("Project.toml")["deps"] |> keys |> y->"GLMakie" in y
    Pkg.rm("GLMakie")
end

Pkg.instantiate()
Pkg.precompile()
