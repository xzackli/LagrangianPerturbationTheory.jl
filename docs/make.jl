using LagrangianPerturbationTheory
using Documenter

DocMeta.setdocmeta!(LagrangianPerturbationTheory, :DocTestSetup, :(using LagrangianPerturbationTheory); recursive=true)

makedocs(;
    modules=[LagrangianPerturbationTheory],
    authors="Zack Li",
    repo="https://github.com/xzackli/LagrangianPerturbationTheory.jl/blob/{commit}{path}#{line}",
    sitename="LagrangianPerturbationTheory.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://xzackli.github.io/LagrangianPerturbationTheory.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/xzackli/LagrangianPerturbationTheory.jl",
    devbranch="main",
)
