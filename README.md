# BuildingGeometry

This package contains Julia code that discretizes the surface of a building.

In wind tunnel tests, the aerodynamic loads are usually obtained from the pressure distribution. Unfotunately, the pressure distribution is obtained from discrete pressure taps located along the surface of the building. To calculate the load, these discrete pressure values should be interpolated. For most cases, except for very simple surface geometries and pressure tap distribution, this is a complex and difficult operation.

This package implements this interpolation by using regions of influence. A region of influence is the zone of the surface closer to a pressure tap than any other tap. What makes this tricky is that the surface might be curved.

The here is to divide the total building surface into surfaces where pressure varies smoothly. For instance a corner separates two faces but a round building could be considered to have a smooth variation of pressure and thus should be treated as a single surface.

Another complication is that the surface has two sides. One of the sides we will call external and it is on the side of the outward normal of the surface. Sometimes only one side has pressure taps. In this case, the internal pressure should be calculated by some other means, using a Wind Code for instance. Other times, the pressure is measured on both sides and the problem is how to combine the external and internal meshes.

The pressure distribution is important but the force distribution is important as well. This package also implements tools that help calculate forces on a surface, slice the surface to obtain loads per slice (usually a floor in tall buildings).


