using ComminWeath
using Documenter

DocMeta.setdocmeta!(ComminWeath, :DocTestSetup, :(using ComminWeath); recursive=true)

makedocs(;
    modules=[ComminWeath],
    authors="Graham Harper Edwards",
    repo="https://github.com/GrahamEdwards/ComminWeath.jl/blob/{commit}{path}#{line}",
    sitename="ComminWeath.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://GrahamEdwards.github.io/ComminWeath.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/GrahamEdwards/ComminWeath.jl",
    devbranch="main",
)
