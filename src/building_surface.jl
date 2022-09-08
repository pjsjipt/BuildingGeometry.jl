
import StaticArrays: SVector
import Meshes: connect, SimpleMesh

struct BuildingSurface{T,Tex,Tin,TMesh}
    tri::Vector{Triangle{3,T,SVector{3,Point{3,T}}}}
    nodes::Vector{NodeInfo{3,T,Tex,Tin,Int}}
    mesh::TMesh
end




function buildsurface(cad::AbstractVector{<:Triangle{3,Float64}},
                      sections::AbstractVector; nointid=-1)

    nsecs = length(sections)

    TriFace = Triangle{3,Float64,SVector{3,Point{3,Float64}}}
    msh = TriFace[]
    nodes = NodeInfo{3,Float64,Int,Int,Int}[]
    # The first section should specify the external side
    # The others, if present are internal sides
    
    # The external and internal sides should span the entire surface.
    # Different internal surfaces shouldn't overlap
    
    # Let's deal with external surface:
    esec = sections[1]
    if !haskey(esec, :points)
        # No external points. Weird...
        evor = [cad[esec.tri]]
        eid = [esec.id[1]]
        eidx = [esec.tri]
    else
        (length(esec.id) != length(esec.points)) &&
            error("`id` should have the same length as `points`")
        evor, eidx = discrsurface(cad, esec.tri, esec.points)
        eid = esec.id
    end
    if nsecs == 1  # No internal surfaces
        for (i, m) in enumerate(evor)
            for (k,t) in enumerate(m)
                An = area(t) .* normal(t)
                push!(msh, Triangle(vertices(t)...))
                push!(nodes, NodeInfo(An, eid[i], nointid, eidx[i][k]))
            end
        end
    else

        Ne = length(evor)
        
        # Let's deal with internal faces
        itri = [collect(sections[i].tri) for i in 2:nsecs]
        etri = esec.tri
        # There should be no overlap
        nisecs = nsecs-1
        for i in 1:nisecs-1
            for k in i+1:nisecs
                if length(itri[i] ∩ itri[k]) > 0
                    error("There should be no overlap in different sections of internal surfaces!")
                end
            end
        end
        itri1 = sort(vcat(itri...))
        etri1 = sort(etri)
        
        if itri1 != etri1
            error("Internal and external surfaces should be the same")
        end
        
        # Let's start building the mesh
        for sec in sections[begin+1:end]
            if !haskey(sec, :points) # No need for voronoi. Just add the triangles
                ivor = [cad[sec.tri]]
                iid = [sec.id[1]]
                iidx = [sec.tri]
            else
                ivor, iidx = discrsurface(cad, sec.tri, sec.points)
                iid = sec.id
            end
            Ni = length(ivor)
            for e in 1:Ne
                for i in 1:Ni
                    intersectmesh!(msh, nodes, eid[e], evor[e], eidx[e],
                                   iid[i], ivor[i], iidx[i])
                end
            end
            
        end
    end
    
    # Lets get the vertices
    ntri = length(msh)
    verts = fill(Point{3,Float64}(0,0,0), ntri * 3)
    conn = fill((0,0,0), ntri)
    for i in 1:ntri
        k = (i-1)*3 + 1
        vtri = vertices(msh[i])
        verts[k] = vtri[1]
        verts[k+1] = vtri[2]
        verts[k+2] = vtri[3]
        conn[i] = (k,k+1,k+2)
    end
    mesh = SimpleMesh(verts, connect.(conn, Triangle))
    
    return BuildingSurface(msh, nodes, mesh)
    
end


