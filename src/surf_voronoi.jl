import GeometryBasics: Point3


project2d(p::SimplePolygon{3,T}, ex, ey) where {T} =
    SimplePolygon([Point2{T}(sum(ex*v), sum(ey*v)) for v in coordinates(p)])
project2d(p::Triangle{3,T}, ex, ey) where {T} =
    Triangle([Point2{T}(sum(ex*v), sum(ey*v)) for v in coordinates(p)]...)
 

function get_basis(n)
    # We want to project the polygons on a surface plane
    # but we still need two two independent vector to form
    # an orthonormal basis. It would be better if the
    # axis are partially aligned with the original
    # 3D reference frame.

    ex = Point3{Float64}(1.0, 0.0, 0.0)
    ey = Point3{Float64}(0.0, 1.0, 0.0)
    ez = Point3{Float64}(0.0, 0.0, 1.0)

    e = (ex, ey, ez)
    dx = sum( ex .* n)
    dy = sum( ey .* n)
    dz = sum( ez .* n)
    # Figure out which axis is more perpendicular to the normal:
    idx = argmin( abs.((dx, dy, dz)) )
    # idx corresponds to the axis more nearly perpendicular.
    # Since in most buildings most faces are vertical or
    # approximately so, the upward direction will be the y axis
    # for the projection.
    
    vy = e[idx] - sum(e[idx] .* n) .* n  # Vector normal to n and e[idx]
    vx = Point3{Float64}(crossprod(vy, n))


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

    ex, ey, ez = get_basis(n⃗)

    # Project each polygon on the surface:
    
end

    

    
