using Documenter, FLPDevelopment

makedocs(;
    modules=[FLPDevelopment],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/DarioSarra/FLPDevelopment.jl/blob/{commit}{path}#L{line}",
    sitename="FLPDevelopment.jl",
    authors="DarioSarra",
    assets=String[],
)

deploydocs(;
    repo="github.com/DarioSarra/FLPDevelopment.jl",
)
