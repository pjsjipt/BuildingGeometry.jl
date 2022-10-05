```@meta
CurrentModule = BuildingGeometry
```

# BuildingGeometry

Documentation for [BuildingGeometry](https://github.com/pjsjipt/BuildingGeometry.jl).

`BuildingGeometry` is a package that handles interpolation of discrete pressure measurements on the surfaces of buildings.

## Introduction

Wind tunnel testing is still the most reliable method of estimating wind loads on structures. Several techniques are used to measure those loads but when dealing buildings, with wide open surfaces, measuring the pressure is one of the best approachs. Not only can the pressure be used to estimate total loads but it can determine how these loads are distributed.

The pressure is measured at discrete points along the surface and the difficulty is interpolating these discrete measurements throughout the surfaces. Except for simple geometries and specific pressure tap distribution this interpolation is not trivial.

Other complications arise such as internal pressure that can be directly measured or estimated by other means (usually Wind Codes or even an indirect measurement).

The approach taken by this package is to assign to each pressure tap a area of influence on each side of the surface. This area of influence is the region of the surface that is closest to one pressure tap in relation to every other tap. The package can handle both flat and curved surfaces. The only restriction is that the pressure is assumed to vary smoothly accross the surface. Thus, when there are corners, the surface should be split into two surfaces

## Algorithm description

Once a surface with associated internal and external pressure taps is available, the package computes the 3D Voronoi diagram of the nodes specified by the pressure taps. This diagram consists of convex polyhedrons that define a volume where everypoint inside this volume is closer to the pressure tap associated with this volume than any other pressure tap on the surface.

Once the Voronoi diagram is available, the surface of the building, represented by a mesh of triangles is split into submeshes, also made up of triangles for each volume. The procedure is repeated with the internal pressure taps and the meshes are intersected.

The result is a mesh of triangles that span the entire surface of the building. Each triangle has an associated external pressure tap and an internal pressure tap, a surface area and an external normal.


## Contents

```@contents
Pages = ["tutorial.md", "docstrings.md"]
```


