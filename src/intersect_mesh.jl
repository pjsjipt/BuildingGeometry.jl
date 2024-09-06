

"""
`intersectmesh(tex, iex, tin, iin)`

Computes the intersection between external and internal triangles that
make up the surface  meshes.
It should be pointed out that this function assumes that both external and internal meshes are overlaid on the same surface mesh.

This function is limited to intersecting one external triangle with one
internal triangle.

# Arguments
 * `msh` List of triangles making up the mesh being built
 * `nodes` `NodeInfo` vector of mesh being built
 * `eid` Identifier of external pressure tap
 * `emsh` Triangles corresponding to an external pressure tap
 * `etri` Index of original mesh triangle that contains `emsh`.
 * `iid` Id of internal pressure tap corresponding to internal triangle
 * `imsh` Triangle of the internal side of the mesh
 * `itri` Index of original mesh triangle that contains `imsh`
 * `etag` Tag of the external side of the surface
 * `itag` Tag of the internal side of the surface
"""
function intersectmesh!(msh::AbstractVector{Tri{3,T}},
                        nodes::AbstractVector{NodeInfo{T,Tuple{Tex,Tin}}},
                        eid::Tex, emsh, etri, iid::Tin,
                        imsh, itri, etag=0, itag=0;
                        atol=atolf(T)) where {T,Tex,Tin}
    
    # This function finds the intersection of two triangular meshes of the same
    # surface (`emsh` and `imsh`) defined by a vector of triangles (`msh`).
    # The two meshes correspond to the two sides of the surface.
    # Not only the meshes of both sides are needed: the index to the original
    # triangle in the cad mesh (`msh`) is also required (`etri` and `itri`).

    # First we need to check if there are any triangles in common between both meshes:
    
    et = sort(unique(etri))  # Independent cad mesh triangles of the external mesh
    it = sort(unique(itri))
    icm = intersect(et, it)

    if length(icm) == 0 # No triangles in common. Do nothing
        return
    end

    for k in icm # We will check each of the common triangles
        for ie in eachindex(etri)
            e = etri[ie]
            if e!= k # Haven't found triangle `k` try next triangle
                continue
            end
            
            te = emsh[ie]
            # Now we will search for the triangle in `itri` that corresponds to `k`
            for ii in eachindex(itri) 
                i = itri[ii]
                if i != k # Not it! Try next
                    continue
                end
                ti = imsh[ii]  
                # `ti` and `te` are part of triangle `msh[k]`.
                # Now we can try to intersect them.
                # Lets compute the intersection
                tri_int = intersect_tri(te, ti, atol=atol) # Let's see if we hava
                        # an intersection
                if length(tri_int) > 0
                    for t in tri_int
                        # Get the centroid of the intersection in m without units
                        tp = centroid(t)
                        # Get the area in m² without units
                        A = area(t)
                        # Let's make shure this is a triangle and not
                        # just an edge:
                        if A > atol^2
                            An = normalarea(t)
                            nn = NodeInfo(An, tp, (eid, iid), (etag, itag))
                            push!(msh, t)
                            push!(nodes, nn)
                        end
                    end
                end
            end
        end 
    end
    return
end


function intersect_tri(tri1::Tri{Dim,T}, tri2::Tri{Dim,T};
                       atol=atolf(T)) where {Dim,T}
    
    # The basic algorithm is to take one triangle as a reference
    # and then take eachside of the triangle as a plane normal
    # to the triangle. This plane will be used to slice the reference
    # triangle.
    v1 = vertices(tri1)
    v2 = vertices(tri2)
    TT = Tri{Dim,T}
        
    n1 = normal(tri1)
    pts = [v2[1], v2[2], v2[3]]


    # First side
    tri1 = let
        pl = Plane(v1[1], (v1[2] - v1[1])×n1)
        cut_with_plane(tri2, pl, atol=atol)
    end
    # Second side.
    # Remember we might have 1 or 2 triangles
    tri2 = let
        pl = Plane(v1[2], (v1[3]-v1[2])×n1)
        tri = TT[]
        for t in tri1
            tx = cut_with_plane(t, pl, atol=atol)
            for t1 in tx
                push!(tri, t1)
            end
        end
        tri
    end

    # Third and last side
    return let
        pl = Plane(v1[3], (v1[1]-v1[3])×n1)
        tri = TT[]
        for t in tri2
            tx = cut_with_plane(t, pl)
            for t1 in tx
                push!(tri, t1)
            end
        end
        tri
    end
    
end
