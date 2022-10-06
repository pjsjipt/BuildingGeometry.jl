
# Tutorial


We will try to use the building in the figure as an example.

![Building](figures/cigarbuilding.svg)

The tutorial will consider two cases:

 1. Simpler building, where only the pressure was measured on the outer surface of the cylinder.
 2. The full case where there are external and internal surfaces with sections where the internal pressure is measured and other sections where it is not.

## A simple cylindrical structure

The building will be cylindrical in shape with a diameter of 30 m and 150 m high. Only external pressure taps exist. They are located along 10 rows equally distributed along the height of the building. Each row is composed of 24 pressure taps totalling 240 pressure taps on the outer surface of the building.

The geometry of th building required by `BuildingGeometry` is a set of faces composed of triangles. Each face is a vector of `Meshes.Triangle`. In this case, the geometry consists of a single cylindrical surface that will be generated in Julia directly

### Defining the geometry of the building

```@example 1
using Meshes
using BuildingGeometry

H = 150.0 # Height
D = 30.0  # Diameter
R = D/2   # Radius

θ = 0.0:15.0:360
nθ = length(θ)
x1 = R * cosd.(θ)
y1 = R * sind.(θ)

p1 = Point.(x1, y1, 0.0)
p2 = Point.(x1, y1, H)

trilst = [Triangle(p1[1], p1[2], p2[2]), Triangle(p1[1], p2[2], p2[1])]

for i in 2:nθ-1
    push!(trilst, Triangle(p1[i], p1[i+1], p2[i+1]))
    push!(trilst, Triangle(p1[i], p2[i+1], p2[i]))
end
```

To plot the meshes, there is the function [`tri2mesh`](@ref) that converts the vector of triangles into a `Meshes.SimpleMesh` that can be visualized with the package [`MeshViz,jl`](https://github.com/JuliaGeometry/MeshViz.jl)

```@example 1
using MeshViz
using Colors
import GLMakie
viz(tri2mesh(trilst));
GLMakie.save("figures/cigarbuilding1.png", GLMakie.current_figure());
```

![Prédio cilindrico](figures/cigarbuilding1.png)


Now, the vector `trilst` contains the only face in the geometry.

### Defining the external pressure taps

We will have 10 equally spaced rows of pressure taps with 24 pressure taps per row.

```@example 1
nr = 10 # Number of rows
dz = H/nr  # Height of each row
zh = range(dz/2, step=dz, length=nr)
θ2 = range(15/2, step=15, length=24)

epts = Point3[]
for z in zh
    for ang in θ2
    	x = R * cosd(ang)
	y = R * sind(ang)
	push!(epts, Point(x, y, z))
    end
end    

viz(epts, color=:red);
viz!(tri2mesh(trilst), color=:gray, alpha=0.3);
GLMakie.save("figures/cigarbuilding2.png", GLMakie.current_figure());
```

![Prédio cilindrico com tomadas](figures/cigarbuilding2.png)

### Discretizing the building

Now that we have the geometry of the building and the position of the pressure tap,
we can discretize the building into regions of influence of each pressure tap. The idea here is to split the surface into triangles with each triangle having two sides. The external side (side 1) and the internal side. Each side will reference the pressure acting on it. In the simplest case, it referes to the index of the pressure tap. Since this is actually a parametric type, it could also be an set of indices of pressure taps and respective weights to compute a weighed average. For now only the simplest case is implemented.

The discretization is done using the function [`buildsurface`](@ref). This function has three arguments:

 1. The geometry of the surface, a vector of triangles
 2. Definition pressure tap positions on different sections of the face. Each section is an element of a vector with a named tuple that has the coordinates of each pressure tap on the section, the index (of the geometry defined in the first argument) that make up the section) and the id of each pressure tap.
 3. The keyword argument specifying the id used for regions that don't have an internal pressure tap.

The function returns a [`BuildingSurface`](@ref) object. This object has an array of triangles that make up the surface, an array of points specifying the coordinates of each triangle (today this is the centroid of each triangle) and an array of [`NodeInfo`](@ref) objects.

The [`NodeInfo`](@ref) objects is the most important field. It stores the influence area of each node times the outer normal. the coordinate of the node and side information as described above. This all that is needed to compute forces and moments. The `BuildingSurface` is just a helper that stores the triangles for a more easy visualization.

In this simple building, only external pressure taps exist. Therefore the internal side (`side==1`) will adopt the `nointid=-1`. There is a single section in this example.

Usually, there are several surfaces that make up the building. After discretizing each surface, they can be merged using method [`mergemeshes`](@ref). 

```@example 1
msh = buildsurface(trilst, # The geometry defined above
                   [(points=epts, # Points defined above
		     tri=1:length(trilst), # Every triangle of the geometry
		     id = 1:240)],
		  nointid = -1); # Use this value for internal pressure tap.

# Now we will try to view the region of influence of each tap
using Colors
cc = distinguishable_colors(240)  # One color for each pressure tap
viz(epts)
ie = nodeside.(msh.nodes, 1)  # Getting the external pressure tap for each triangle
viz!(tri2mesh(msh.tri), color=cc[ie])
GLMakie.save("figures/cigarbuilding3.png", GLMakie.current_figure());
```

![Influence regions of each pressure tap](figures/cigarbuilding3.png)


### Slicing the building

An important result of wind tunnel testing is obtaining a distribution of forces. For a tall building this means forces and moments on each floor. With this in mind, the total mesh of the building must, often, be sliced so that the forces on each floor can be computed. The generic function [`buildingslice`](@ref) has several methods to slice the building. In our simple cylindrical building, we will assume that each floor is 3 m high.

```@example 1
zslices = 0.0:3.0:H  # Boundaries of each slice

slices = buildingslice(msh, zslices);

# Let's try to plot every other floor
viz(epts, color=:black, size=3) # Pressure taps
for i in firstindex(slices):2:lastindex(slices)
    ie1 = nodeside.(slices[i].nodes, 1)
    viz!(tri2mesh(slices[i].tri), color=cc[ie1])
end
GLMakie.save("figures/cigarbuilding4.png", GLMakie.current_figure());
```

![Every other slice](figures/cigarbuilding4.png)


### Computing the forces

With the building discretized into triangles and with each triangle having information on the outside and inside pressure, forces can be easily calculated with the [`NodeInfo`](@ref) obtained with [`buildsurface`](@ref) and [`buildingslice`](@ref).

Given a pressure distribution, the force acting on one side of the surface is given by

$\vec{F} = -\int_S p \:d\vec{A}$

The moment is calculated by

$\vec{M} = -\int_S p \vec{r}\times d\vec{A}$

Using the discretization obtained above,

$\vec{F} = -\sum_{i=1}^{N_{taps}} p_i \cdot \vec{A}_i$

$\vec{M} = -\sum_{i=1}^{N_{taps}} p_i \cdot \vec{r_i} \times \vec{A}_i$

Notice that given a vector with every pressure measurement, the operations above can be represented as a matrix multiplication.

$\left\{\begin{matrix}F_x\\F_y\\F_z\\M_x\\M_y\\_Mz\end{matrix}\right\} = \left[ F_{matrix} \right] \cdot \left\{\begin{matrix}p_1\\p_2\\p_3\\\vdots\\p_{N_{taps}}\end{matrix}\right\}$

The $\left[F_{matrix}\right] matrix is sparse. The number of rows is the number of triangles (or nodes in more general cases) and the number of columns is the number of pressure taps. To compute this matrix for a discretized surface, use the [`addforcecontrib!`](@ref) method or [`forcematrix`](@ref). The `forcematrix` method allocates memory for the matrix and calls `addforcecontrib!` method. The methods are defined in this way because there might be contributions from both sides of the surface and each contribution should be added independently to result in the full force matrix.

```@example 1
# Remember we have 240 pressure taps!
Fbase = forcematrix(240, msh, (1,2,3,4,5,6); sgn=1, side=1, point=Point(0,0,0))
println("Dimensions of `Fbase`: $(size(Fbase))")
```

The first parameter is the number of columns. This is actually the number of pressure taps  used. The second parameter is the mesh. The third parameter specifies which components of the force should be calculated:

 1. $F_x$
 2. $F_y$
 3. $F_z$
 4. $M_x$
 5. $M_y$
 6. $M_z$

The keyword argument `sgn` multiplies each matrix element. In the case of internal pressure, usually this argument is -1. But it can also account for scaling factors or different reference wind velocity. The `side` keyword argument specifies which side of the face is contributing to a force. In a more general case, `forcematrix` would be called sith `sgn=1` and `side=1` then `addforcecontrib!` would be called with `sgn=-1` and `side=2`. Since we only have external pressure taps, in the example above, `sgn=1` and `side=1`.

The moments of each triangle is calculated with respect to the point specified by the `point` keyword argument.

The matrix above can be used to calculate the loads on the foundation of the building.


### Calculating the forces on each slice

Usually the structure designer wants a load distribution. For instance on a tall building, the finite element software usually requires the loads on each floor. There are other methods for [`addforcecontrib!`](@ref) and [`forcematrix`](@ref) for dealing with multiple meshes simultaneously.

In this case, the force matrix, when applied to the pressure measurements calculates the forces on each mesh. The `interleaved` keyword argument specifies how the forces are numbered in sequence:

 * `interleaved=false`: the each component of the force is numbered in sequence. For example if the argument `forces=(1,2,6)` the order of the forces is $F_{x,1}$, $F_{x,2}$, $\ldots$, $F_{x,N}$, $F_{y,1}$, $F_{y,2}$, $\ldots$, $F_{y,N}$, $M_{z,1}$, $M_{z,1}$, $\ldots$, $M_{z,N}$ where $N$ is the number of floors.
 * `interleaved=false`: The forces are numbered per floor (or mesh), $F_{x,1}$, $F_{y,1}$, $M_{z,1}$, $F_{x,2}$, $F_{y,2}$, $M_{z,2}$, $\ldots$, $F_{x,N}$, $F_{y,N}$, $M_{z,N}$.


```@example 1

Fslices = forcematrix(240, slices, (1,2,6); sgn=1, side=1, point=Point(0,0,0))
println("Dimensions of `Fslices`: $(size(Fslices))")
```


## A more complex building


This new building has all the percs shown in the original figure. It has an internal division dividing the cylinder in two halves. One halve is isolated and therefore there are not internal pressure taps (as was the case in the simple building above). But the other half has both external pressure taps and internal pressure taps.


This model has two surfaces:
 1. The cylindrical surface with external pressure taps all over (the same as the simple building above). This surface has two sections:
    * One section (half of the cylinder) has both internal and external pressure taps.
    * The second section, the other half has only external pressure taps
 2. The flat surface in that splits the cylinder in half with pressure taps on one side only (we will call this side the exterior).

### Defining the geometry

We will start out with the original geometry of the cicrcular building.

```@example 2
using Meshes, MeshViz
import GLMakie
using BuildingGeometry

H = 150.0 # Height
D = 30.0  # Diameter
R = D/2   # Radius

θ = 0.0:15.0:360
nθ = length(θ)
x1 = R * cosd.(θ)
y1 = R * sind.(θ)

p1 = Point.(x1, y1, 0.0)
p2 = Point.(x1, y1, H)

face1 = [Triangle(p1[1], p1[2], p2[2]), Triangle(p1[1], p2[2], p2[1])]

for i in 2:nθ-1
    push!(face1, Triangle(p1[i], p1[i+1], p2[i+1]))
    push!(face1, Triangle(p1[i], p2[i+1], p2[i]))
end


pf1 = Point(R, 0, 0)
pf2 = Point(-R, 0, 0)
pf3 = Point(-R, 0, H)
pf4 = Point(R, 0, H)
face2 = [Triangle(pf1, pf2, pf3), Triangle(pf1, pf3, pf4)]

viz(tri2mesh(face1), color=:blue);
viz!(tri2mesh(face2), color=:red);
GLMakie.save("figures/building1.png", GLMakie.current_figure());
```

![Prédio](figures/building1.png)

### Defining the external and internal pressure taps

The external pressure taps are the same as the simple building. But now we need to add internal pressure taps and external pressure taps on the second face.


```@example 2
nr = 10 # Number of rows
dz = H/nr  # Height of each row
zh = range(dz/2, step=dz, length=nr)
θ2 = range(15/2, step=15, length=24)

epts1 = Point3[]

for z in zh
    for ang in θ2
    	x = R * cosd(ang)
	y = R * sind(ang)
	push!(epts1, Point(x, y, z))
    end
end    

# Internal nodes of face 1
nri = 3
dzi = H / nri
zhi = range(dzi/2, step=dzi, length=nri)
θi = range(15.0, step=30, length=6)
ipts1  = Point3[]

for z in zhi
    for ang in θi
    	x = R * cosd(ang)
	y = R * sind(ang)
	push!(ipts1, Point(x, y, z))
    end
end



# External nodes of face 2

nx2 = 3
dx2 = D/nx2
x2 = range(-R+dx2/2, step=dx2, length=nx2)
epts2 = Point3[]

for z in zhi
    for x in x2
    	push!(epts2, Point(x, 0.0, z))
    end
end




viz(epts1,color=:red)
viz!(ipts1, color=:blue)
viz!(epts2, color=:green)

viz!(tri2mesh(face1), color=:gray, alpha=0.3);
viz!(tri2mesh(face2), color=:gray, alpha=0.3)

GLMakie.save("figures/building2.png", GLMakie.current_figure());
```

![Building with pressure taps](figures/building2.png)
