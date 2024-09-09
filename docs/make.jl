using Documenter
using Colors, GLMakie
using BuildingGeometry

makedocs(
    sitename = "BuildingGeometry",
    format = Documenter.HTML(prettyurls=false),
    modules = [BuildingGeometry]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
