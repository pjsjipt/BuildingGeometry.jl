


"""
`intersectmesh(tex, iex, tin, iin)`

Computes the intersection between external and internal triangles that
make up the surface  meshes.
It should be pointed out that this function assumes that both external and internal meshes are overlaid on the same surface mesh.

This function is limited to intersecting one external triangle with one
internal triangle.

# Arguments
 * `msh` List of triangles making up the mesh being built
 * `nodes` `NodeInfo` vector of mesh being build
 * `eid` Identifier of external pressure tap
 * `emsh` Triangles corresponding to an external pressure tap
 * `etri` Index of original mesh triangle that contains `emsh`.
 * `iid` Id of internal pressure tap corresponding to internal triangle
 * `imsh` Triangle of the internal side of the mesh
 * `itri` Index of original mesh triangle that contains `imsh`
 * `etag` Tag of the external side of the surface
 * `itag` Tag of the internal side of the surface
"""
function intersectmesh!(msh::AbstractVector{Tri},
                        nodes::AbstractVector{NodeInfo{3,T,Tuple{Tex,Tin}}},
                        eid::Tex, emsh, etri, iid::Tin,
                        imsh, itri, etag=0, itag=0; rtol=1e-8) where {Tri,T,Tex,Tin}
    et = sort(unique(etri))  # Independent cad mesh triangles of the external mesh
    it = sort(unique(itri))
    icm = intersect(et, it)

    if length(icm) == 0
        return
    end
    lidx = lastindex(msh)
    for k in icm
        for ie in eachindex(etri)
            e = etri[ie]
            if e!= k
                continue
            end
            
            te = emsh[ie]
            for ii in eachindex(itri)
                i = itri[ii]
                if i != k
                    continue
                end
                ti = imsh[ii]
                # Lets compute the intersection
                pts = intersect_tri(te, ti; rtol=rtol)
                if length(pts) > 2
                    for l in firstindex(pts)+1:lastindex(pts)-1
                        new_tri = Triangle(pts[begin], pts[l], pts[l+1])
                        tp = centroid(new_tri)
                        A = area(new_tri)
                        if A > 0
                            An = A .* normal(new_tri)
                            nn = NodeInfo(An, tp, (eid, iid), (etag, itag))
                            push!(msh, new_tri)
                            push!(nodes, nn)
                        end
                    end
                end
                
                
            end
        end
    end
    return
end


function intersect_tri(tri1, tri2; rtol=1e-8)
    v1 = vertices(tri1)
    v2 = vertices(tri2)

    A1 = area(tri1)
    A2 = area(tri2)
    Lref = sqrt(max(A1,A2))
    atol = rtol*Lref
    

    n1 = normal(tri1)
    pts = [v2[1], v2[2], v2[3]]

    pts = cut_with_plane(pts, v1[1], (v1[2] - v1[1]) × n1; atol=atol)
    if length(pts) == 0 && 
        return pts
    end
    pts = cut_with_plane(pts, v1[2], (v1[3] - v1[2]) × n1; atol=atol)
    if length(pts) == 0 && 
        return pts
    end
    pts = cut_with_plane(pts, v1[3], (v1[1] - v1[3]) × n1; atol=atol)
    return pts
    
end
