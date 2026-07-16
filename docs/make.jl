using BoreholeResistance
using Documenter

# Make `using BoreholeResistance` available to every doctest in docstrings and pages.
DocMeta.setdocmeta!(
    BoreholeResistance,
    :DocTestSetup,
    :(using BoreholeResistance);
    recursive = true,
)

makedocs(;
    modules = [BoreholeResistance],
    authors = "Gabriel-Dion <dion.gabriel100@gmail.com>",
    sitename = "BoreholeResistance.jl",
    format = Documenter.HTML(;
        canonical = "https://GHE-jl.github.io/BoreholeResistance.jl",
        edit_link = "master",
        assets = String[],
        mathengine = Documenter.KaTeX(),
        sidebar_sitename = false,
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Modeling theory" => [
            "Resistance network" => "theory/overview.md",
            "Fluid convective resistance" => "theory/fluid.md",
            "Pipe conductive resistance" => "theory/pipe.md",
            "Borehole (grout) resistance" => "theory/borehole.md",
            "Effective resistance" => "theory/effective.md",
        ],
        "Water properties" => "properties.md",
        "API reference" => "api.md",
        "References" => "references.md",
    ],
    # Keep the build strict so broken cross-references or missing docstrings fail CI.
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/GHE-jl/BoreholeResistance.jl",
    devbranch = "master",
)
