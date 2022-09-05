

let

    pts = Point.([(-0.5,-0.5,0),(0.5,-0.5,0),(-0.5,0.5,0),(0.5,0.5,0)])
    tri = [Triangle((-1,-1,0), (1,-1,0), (1,1,0)), Triangle((-1,-1,0),(1,1,0),(0,1,0))]
    vor = discrsurface(tri, 1:2, pts)
    
           
end
