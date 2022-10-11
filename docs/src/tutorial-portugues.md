
# Tutorial traduzido para o Português

Vamos tentar trabalhar com o modelo do edifício abaixo como exemplo

![Building](figures/cigarbuilding.svg)

Este tutorial é dividido em duas partes:

 1. Edifício simples onde apenas a pressão externa do edifício cilíndrivo é medida.
 2. O caso completo da figura acima onde existem 2 fases com seções onde a pressão interna é medida ou não.


## Um edifício cilíndrico simples

O edifício é um cilindro com diâmetro de 30 m e 150 m de altura. Apenas a pressão externa é medida. As tomadas de pressão estão localizadas em 10 linhas igualmente distribuídas ao longo da altura do edifício. Cada linha é composta de 24 tomadas de pressão, totalizando 240 tomadas de pressão externa em uma única face.

O pacote `BuildingGeometry` trabalha com superfícies composta por triângulos. Cada superfície é um vetor de `Meshe.Triangle`. Neste caso, a geometria se reduz a uma única face cilíndrica que será gerada em Julia diretamente.

### Definição da geomtetria do edifício

```@example 3
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

Para se visualizar a geometria, existe a função [`tri2mesh`](@ref) que converte o vetor de triângulos em um objeto `Meshes.SimpleMesh` que pode ser visualizado com o pacote [`MeshViz,jl`](https://github.com/JuliaGeometry/MeshViz.jl).

```@example 3
using MeshViz
using Colors
import GLMakie
viz(tri2mesh(trilst));
GLMakie.save("figures/cigarbuilding1.png", GLMakie.current_figure());
```

![Prédio cilindrico](figures/cigarbuilding1.png)

Agora, o vetor `trilst` contém a única face da geometria.

### Definindo as tomadas de pressão externals

Estas estarão distribuídas ao longo de 10 linhas igualmente espaçadas com 24 tomadas de pressão por linha.

```@example 3
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

### Discretização do edifício

Agora que temos a geometria do edifício e as posições das tomadas de pressão, podemos discretizar a superfície do edifício em zonas de influência de cada tomada de pressão. A idéia aqui é dividir a superfície em triângulos onde cada triângulo tem dois lados. O externo (lado 1) e o interno (lado 2). Cada lado terá informação sobre a pressão agindo. Em geral isso pode ser uma média ponderada das pressões medidas mas no caso mas simples é a tomada de pressão mais próxima. Neste exemplo, apenas a pressão externa é usada e cada triângulo aponta para uma única tomada de pressão.

A discretização é feita usando a função [`buildsurface`](@ref). Esta função tem dois argumentos:

 1. A geometria da superfície
 2. Definição da posição das tomadas de pressão em cada face. Cada seção é um elemento de um vetor com as coordenadas da tomada de pressão (campo `points`), 


