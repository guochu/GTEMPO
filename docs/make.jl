using Documenter
using GTEMPO

DocMeta.setdocmeta!(GTEMPO, :DocTestSetup, :(using GTEMPO); recursive=true)

makedocs(;
    modules=[GTEMPO],
    authors="Guo Chu, Zhijie Sun",
    sitename="GTEMPO.jl",
    warnonly=[:missing_docs],
    format=Documenter.HTML(;
        canonical="https://example.com/GTEMPO.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Practice guide" => "practice.md",
        "Grassmann lattices" => "grassmann_lattice.md",
        "Tensor and Grassmann conventions" => "tensor_and_grassmann_conventions.md",
        "Internals" => "internals.md",
        "Tutorials" => "tutorials.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/example/GTEMPO.jl",
)
