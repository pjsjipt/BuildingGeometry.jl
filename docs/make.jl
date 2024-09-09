using GLMakie
using BuildingGeometry
using Documenter

DocMeta.setdocmeta!(BuildingGeometry, :DocTestSetup, :(using BuildingGeometry); recursive=true)

makedocs(;
    modules=[BuildingGeometry],
    authors="Paulo Jabardo <pjabardo@ipt.br>",
    repo="https://github.com/pjsjipt/BuildingGeometry.jl/blob/{commit}{path}#{line}",
    sitename="BuildingGeometry.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
