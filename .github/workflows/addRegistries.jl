#= The MIT License (MIT)

Copyright (c) 2018-2020 GitHub, Inc., David Anthoff and contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
=#
using Pkg

function general_registry_location()
    general_registry_dir = joinpath(DEPOT_PATH[1], "registries", "General")
    registry_toml_file = joinpath(general_registry_dir, "Registry.toml")
    return general_registry_dir, registry_toml_file
end

function general_registry_exists()
    general_registry_dir, registry_toml_file = general_registry_location()
    if !isdir(general_registry_dir)
        return false
    elseif !isfile(registry_toml_file)
        return false
    else
        return true
    end
end

function add_general_registry()
    @info("Attempting to clone the General registry")
    general_registry_dir, registry_toml_file = general_registry_location()
    rm(general_registry_dir; force = true, recursive = true)
    Pkg.Registry.add("General")
    isfile(registry_toml_file) || throw(ErrorException("the Registry.toml file does not exist"))
    return nothing
end

function main(; n = 10, max_delay = 120)
    VERSION >= v"1.5-" || return
    print("Running local add General")
    if general_registry_exists()
        @info("The General registry already exists locally")
        return
    end

    delays = ExponentialBackOff(; n = n, max_delay = max_delay)
    try
        retry(add_general_registry; delays = delays)()
        @info("Successfully added the General registry")
    catch ex
        msg = "I was unable to add the General registry. However, the build will continue."
        @error(msg, exception=(ex,catch_backtrace()))
    end

    return
end

main()
