using GLMakie
using BuildingGeometry
using Documenter

DocMeta.setdocmeta!(BuildingGeometry, :DocTestSetup, :(using BuildingGeometry); recursive=true)

makedocs(;
         modules=[BuildingGeometry],
         authors="Paulo Jabardo <pjabardo@ipt.br>",
         sitename="BuildingGeometry.jl",
         format=Documenter.HTML(prettyurls = haskey(ENV, "CI")),
         pagesonly=true,
         draft=false,
         pages=["Home" => "index.md"])


#         repo="https://github.com/pjsjipt/BuildingGeometry.jl/blob/{commit}{path}#{line}",
