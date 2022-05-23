#import GeometryBasics: Point, Point3, Triangle
import Meshes: Point, Point3, Point2
using LinearAlgebra
import StaticArrays
import GeneralPolygonClipper as gpc
using VoronoiCells


struct SProject2d end


project2d(p::ConvexPolygon{3,T}, ex, ey) where {T} =
    ConvexPolygon([Point2{T}(sum(ex*v), sum(ey*v)) for v in coordinates(p)])
project2d(p::Triangle{3,T}, ex, ey) where {T} =
    Triangle([Point2{T}(sum(ex*v), sum(ey*v)) for v in coordinates(p)]...)

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

normal1(x) = (x[2]-x[1]) × (x[3]-x[1])

    


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

    # Calculate the mean normal
    nrm = Vec3{Float64}(0,0,0)
    for f in faces
        for tri in f
            n = normal1(tri)
            nrm += A*n
        end
    end

    nrm = nrm ./ hypot(nrm...)
    
    # Project points and faces on a plane normal to `nrm`
    e = get_basis(nrm) # Basis


    P₀ = Point3(0,0,0)
    # Project the triangles
    faces2d = [[project_face(e, P₀, tri3) for tri3 in ff] for ff in faces]

    # Project the pressure taps
    pe2d = [Point2( (p - P₀) ⋅ e[1], (p-P₀)⋅e[2] ) for p in pext]
    pi2d = [Point2( (p - P₀) ⋅ e[1], (p-P₀)⋅e[2] ) for p in pint]

    # of the points
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

    # Boundaing box
    rect = Rectangle(Point2(xmin, ymin), Point2(xmax, ymax))

    # Voronoi tesselation of external points
    tessex = voronoicells(pe2d, rect)

    # Let's go through every Voronoi cell and intersect it
    for i in eachindex(pe2d)
        vpts = tessex.Cells[i] # Voronoi cell
        (xmin1,ymin1), (xmax1,ymax1) = boundingbox(vpts) 
        for face in faces2d
            for tri in face
                (xmin2, ymin2), (xmax2, ymax2) = boundingbox2(tri)
                xposs = xmin1 > xmax2 || xmin2 
            end
        end
        
    
    
end

                      
function project_face(e, p, tri3)
    v1 = Point2( (tri3[1]-p) ⋅ e[1], (tri3[1]-p) ⋅ e[2])
    v2 = Point2( (tri3[2]-p) ⋅ e[1], (tri3[2]-p) ⋅ e[2])
    v3 = Point2( (tri3[3]-p) ⋅ e[1], (tri3[3]-p) ⋅ e[2])
    return TriangleFace(v1, v2, v3)
    
end


    
