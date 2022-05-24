#import GeometryBasics: Point, Point3, Triangle
import GeometryBasics: Point, Point3, Point2, Vec, Vec2, Vec3, TriangleFace, Rectangle
using LinearAlgebra
import StaticArrays
import GeneralPolygonClipper as gpc
using VoronoiCells


struct SProject2d end



function boundingbox2(pts)
    i = firstindex(pts)
    xmin = pts[i][1]
    ymin = pts[i][2]

    xmax = xmin
    ymax = ymin

    for p in pts
        x,y = p
        if x < xmin
            xmin = x
        elseif x > xmax
            xmax = x
        end
        if y < ymin
            ymin = y
        elseif y > ymax
            ymax = y
        end
    end
    return (xmin, ymin), (xmax, ymax)
    
end

function boundingbox3(pts)
    
    i = firstindex(pts)
    xmin = pts[i][1]
    ymin = pts[i][2]
    zmin = pts[i][3]

    xmax = xmin
    ymax = ymin
    zmax = zmin
    
    for p in pts
        x,y,z = p
        if x < xmin
            xmin = x
        elseif x > xmax
            xmax = x
        end
        if y < ymin
            ymin = y
        elseif y > ymax
            ymax = y
        end
        if z < zmin
            zmin = z
        elseif z > zmax
            zmax = z
        end
        
    end
    return (xmin, ymin, zmin), (xmax, ymax, zmax)
    
end


function get_basis(n::Vec{3,T}) where {T}
    # We want to project the polygons on a surface plane
    # but we still need two two independent vector to form
    # an orthonormal basis. It would be better if the
    # axis are partially aligned with the original
    # 3D reference frame.

    ex = Vec{3,T}(1, 0, 0)
    ey = Vec{3,T}(0, 1, 0)
    ez = Vec{3,T}(0, 0, 1)

    e = (ex, ey, ez)
    dx = ex ⋅ n 
    dy = ey ⋅ n
    dz = ez ⋅ n 

    # If the surface is almost horizontal idxmax=3
    idxmax = argmax( abs.((dx, dy, dz)) )
    if idxmax==3
        idx = 2 # Y axis will be vertical
    else
        # Figure out which axis is more perpendicular to the normal:
        idx = argmin( abs.((dx, dy, dz)) )
    end
    # idx corresponds to the axis more nearly perpendicular.
    # Since in most buildings most faces are vertical or
    # approximately so, the upward direction will be the y axis
    # for the projection.
    
    vy = e[idx] - (e[idx] ⋅ n) * n  # Vector normal to n and e[idx]
    vx = vy × n


    # return the new *orthonormal* basis
    return vx ./ hypot(vx...), vy ./ hypot(vy...), n ./ hypot(n...)
    
end
    

""" 
`mean_mesh_normal(msh)`

Compute the mean normal of a mesh
"""
function mean_mesh_normal(msh)

    np = length(msh) # Number of polygons (or triangles...)

    nrm <- sum(normal(p) for p in msh)
    n⃗ = nrm ./ hypot(nrm...) # Unit vector

    return get_basis(n⃗)
    
end

normal1(x) = (x[2]-x[1]) × (x[3]-x[1]) / 2

    


"""
`discrsurface(type, faces, pext, fext, pint, fint)`

Building surface discretization. 

This function discretizes a face of a building around 
 internal and external pressure taps. A face is made up of subfaces
which are made up of triangles.

 * `type` A type specifying the algorithm used
 * `faces` a collection of faces
 * `pext` External pressure taps
 * `fext` Index of faces where the pressure taps are located
 * `pint` Internal pressure taps
 * `fint` Index of faces where internal taps are located. Should be a subset of `fext`
"""
function discrsurface(::Type{SProject2d}, faces, pext, fext,  pint, fint)
    # Check if `fint` is an subset of `fext`
    for f in fint
        if f ∉ fext
            error("External surfaces should include every internal surface!")
        end
    end

    # Calculate the mean normal of all external faces
    nrm = Vec3{Float64}(0,0,0)
    for f in faces[fext]
        for tri in f
            n = normal1(tri)
            nrm += n
        end
    end

    nrm = nrm / norm(n)
    
    # Project points and faces on a plane normal to `nrm`
    e = get_basis(nrm) # Basis


    P₀ = Point3(0,0,0)

    # Discretize external surface:
    tri_ext = partition_faces(pext, faces, fext, e, P₀)

    # Discretize internal surfaces, if there are any
    if length(pint) > 0
        # Get the mesh
        tri_int = partition_faces(pint, faces, fint, e, P₀)

        # Now we should intersec
    end
    
end

"""
`partition_faces(pts, facelist, fidx, e, p0)`

Partition a collection of faces in regions corresponding to a collection of points.

The faces are collections of triangles and they are projected into a plane defined
by point `p0` and normal `e[3]`. This plane has two direction vectors:

 * `e[1]` that correspnds to the 2d axis direction vector î
 * `e[2]` that correspnds to the 2d axis direction vector ĵ

Once the points and faces are projected into the plane, a Voronoi tesselation
is carried out with respect to the points and each Voronoi cell is intersected
with each triangle that compose the faces. This intersection is represented
as a set of triangles.

## Arguments

 * `pts` Vector of points
 * `facelist` A vector containing every face in the model
 * `fidx` Vector with indices of faces that should be processed
 * `e` Basis defining the projection plane
 * `p0` Origin of the 2d coordinate system on the projection plane

## Return value

This function returns a vector with a tuple containing a triangle, the index of
the point, the index of the face that this triangle belongs to and the index of the 
model triangle to which the triangle belongs.

"""
function partition_faces(pts, facelist, fidx, e, p0)

    # Let's work with the necessary faces only!
    faces = facelist[fidx]

    # A plane is specified by a point and a normal. Just use the origin as the point
    P₀ = p0
    # Project the triangles
    faces2d = [[project_triangle(e, P₀, tri3) for tri3 in ff] for ff in faces]

    # Project the points on thge same plane
    pts2 = [Point2( (p - P₀) ⋅ e[1], (p-P₀)⋅e[2] ) for p in pts]

    # Get a bounding box of all the projected surfaces to create the Voronoi tesselation
    xmin, ymin = faces2d[1][1][1]
    xmax, ymax = xmin, ymin
    for f in faces2d
        for tri in f
            for i in 1:3
                x,y = tri[i]
                if x < xmin
                    xmin = x
                elseif x > xmax
                    xmax = x
                end
                if y < ymin
                    ymin = y
                elseif y > ymax
                    ymax = y
                end
            end
        end
    end
    Δx = xmax - xmin; Δy = ymax - ymin
    
    # Bounding box. Let's make it bigger to ensure everything fits...
    rect = Rectangle(Point2(xmin-Δx, ymin-Δy), Point2(xmax+Δx, ymax+Δy))

    # Voronoi tesselation of external points
    tessex = voronoicells(pts2, rect)

    # Array to store the 2d triangles, the face index and triangle index.
    tri2d = Tuple{TriangleFace{Point2{Float64}}, Int, Int, Int}[]
    
    # Let's go through every Voronoi cell and intersect it with every triangle
    # in the faces
    for i in eachindex(pts)
        vpts = tessex.Cells[i] # Voronoi cell
        (xmin1,ymin1), (xmax1,ymax1) = boundingbox2(vpts)
        # Create a GenericPolygonClipper.jl polygon
        pgpc = gpc.GPCPolygon([false], [[gpc.Vertex(v[1], v[2]) for v in vpts]])
        for (k, face)   in enumerate(faces2d)
            for (j, tri) in enumerate(face)
                (xmin2, ymin2), (xmax2, ymax2) = boundingbox2(tri)
                if !(xmax2 < xmin1 || xmax1 < xmin2 ||
                    ymax2 < ymin1 || ymax1 < ymin2)
                    # There *might* be an intersection!
                    # If there is an intersection, we want it decomposed in triangle
                    # strips.
                    trigpc = gpc.GPCPolygon([false],
                                            [[gpc.Vertex(tri[i][1], tri[i][2]) for i in 1:3]])
                    # Now we carry out the intersection between the Voronoi cell
                    # and the face triangle
                    trilst = gpc.intersect_strip(pgpc, trigpc)
                    for s in trilst
                        for l in 1:length(s)
                            vv = s[l]
                            tri_new = TriangleFace(Point2{Float64}(vv[1].x, vv[1].y),
                                                   Point2{Float64}(vv[2].x, vv[2].y),
                                                   Point2{Float64}(vv[3].x, vv[3].y))
                            push!(tri2d, (tri_new, i, fidx[k], j))
                        end
                    end
                    
                    
                end
                
            end
        end
    end
    # Project back to 3d.
    # And remember that the orientation might have been inverted
    return [ (project_back(e, P₀, x[1], facelist[x[3]][x[4]]),
              x[2], x[3], x[4]) for x in tri2d ]
    
    
    
  

end


"""
`project_triangle(e, p, tri3)`

Project a 3d triangle into a plane defined by basis `e` with origin `p`.

## Arguments

 * `e` Tuple containing the basis. `e[1]` and `e[2]` are the vectors definind the plane and `e[3]` is the vector normal to the plane
 * `p` 3d point which is used as origin in the plane
 * `tri3` 3d Triangle that should be projected.

## Return value

A `2d TriangleFace` object containing the projected triangle.
"""                      
function project_triangle(e, p, tri3)
    v1 = Point2( (tri3[1]-p) ⋅ e[1], (tri3[1]-p) ⋅ e[2])
    v2 = Point2( (tri3[2]-p) ⋅ e[1], (tri3[2]-p) ⋅ e[2])
    v3 = Point2( (tri3[3]-p) ⋅ e[1], (tri3[3]-p) ⋅ e[2])
    return TriangleFace(v1, v2, v3)
    
end


    
"""
`project_back(e, p, tri2d, tri0)`

Project back a 2d triangle into the 3d surface of the model.


This function projects a 2d triangle back into 3d space where the triangle 
belongs to the surface defined by 3d triangle `tri0`. 

The 2d plane where `tri2d` is defined has as origin the point `p` (3d) and
unit direction vectors `î = e[1]` and `ĵ = e[2]`.

The function will reorder the output vertices so that it has the same orientation
as the 3d triangle `tri0`.

## Arguments
 * `e` Tuple with direction vectors defining the projection plane
 * `p` 3d point which is used as origin in the 2d projection plane
 * `tri2d` The 2d triangle that should be projected back
 * `tri0` triangle that defines the 3d plane where the output triangle should be located.

## Return value

A `TriangleFace` corresponding to `tri2d` in the plane defined by `tri0`.
"""
function project_back(e, p, tri2d, tri0)
    # Any point P in the same plane of a triangle
    # can be computed as P - P₁ = α⋅u + β⋅v
    # Where:
    #  * P₁ is the first vertex of the triangle
    #  * u is the vector from P₁ to P₂
    #  * v is the vector from P₁ to P₃
    # In the projected plane (basis `e`), this results
    # in a system of equations M = [ux vx; uy vy]
    # M⋅[α; β] = [x-x₁, y-y₁]
    # Inverting M we can calculate α and β and now
    # we can use the original definintion above to compute P.

    # Project first vertex on the plane
    vx₁ = (tri0[1]-p) ⋅ e[1]
    vy₁ = (tri0[1]-p) ⋅ e[2]

    u = tri0[2] - tri0[1]
    v = tri0[3] - tri0[1]
    ux = u ⋅ e[1]
    uy = u ⋅ e[2]

    vx = v ⋅ e[1]
    vy = v ⋅ e[2]
    
    Δ = ux * vy - uy * vx

    A11 =  vy / Δ; A12 = -vx / Δ
    A21 = -uy / Δ; A22 =  ux / Δ

    δx = tri2d[1][1] - vx₁
    δy = tri2d[1][2] - vy₁
    α = A11*δx + A12*δy
    β = A21*δx + A22*δy
    p1 = tri0[1] + α * u + β * v

    δx = tri2d[2][1] - vx₁
    δy = tri2d[2][2] - vy₁
    α = A11*δx + A12*δy
    β = A21*δx + A22*δy
    p2 = tri0[1] + α * u + β * v

    δx = tri2d[3][1] - vx₁
    δy = tri2d[3][2] - vy₁
    α = A11*δx + A12*δy
    β = A21*δx + A22*δy
    p3 = tri0[1] + α * u + β * v

    # Check orientation of the nodes.
    # If the normals are in the same direction, keep the order p1, p2, p3.
    # Otherwise invert: p1, p3, p2.
    if (u × v) ⋅ ( (p2-p1)×(p3-p1) ) > 0  
        return TriangleFace(p1, p2, p3)
    else
        return TriangleFace(p1, p3, p2)
    end
    
        
    
    
end
