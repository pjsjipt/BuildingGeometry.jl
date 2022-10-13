
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

 1. A geometria da superfície (argumento `cad`)
 2. Definição da posição das tomadas de pressão em cada seção da face. A primeir seção é necessariamente a a face externa da superfície. Cada seção é um elemento de um vetor com as coordenadas da tomada de pressão (campo `points`), os índices dos triângulos que fazem formam a seção (campo `tri`), a identificação da tomada de pressão (campo `id`) e um tag associado a esta seção (campo `tag`).

Caso o campo `id` ou `tag` não esteja especificado em alguma seção, os valores padrão especificados pelos argumentos palavra chave (*keyword arguments*) `nointid` e `tag` são utilizados.


A função returna um objeto [`BuildingSurface`](@ref). Este objeto tem um vetor de triângulos que formam a superfície, um vetor de pontos com o centróide de cada triângulo e um vetor com objetos [`NodeInfo`](@ref).

Estes objetos [`NodeInfo`](@ref) representam o campo mais importante. Este campo armazena a área de influência vezes a normal externa de cada nó, a coordenada do nó, informação sobre o que está em cada lado da superfície e um tag para identificar onde está o nó. Isto contém tudo para calcular forças e momentos. A estrutura `BuildingSurface` simplesmente auxilia na plotagem da superfície e resultados.

Neste exemplo simples, existem apenas tomadas de pressão externas e uma única superfície. O lado interno (`side == 2`) terá como pressão interna o valor -1 (especificado via argumento palavra chave `nointid == -1` ou diretamente pelo campo `id` descrito acima.

Em problemas mais complexos (e realistas...), mais de uma superfície existe. Elas podem ser juntadas chamando o método [`mergemeshes`](@ref).


```@example 3
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


### Fatiando o edifício

Um importante resultado do teste em túnel de vento é obter a distribuição de forças. No caso de edifícios altos isso significa forças e momentos agindo em cada pavimento. Tendo isso em mente, a malha total do edifício, frequentemente deve ser fatiado para que a cada pavimento exista uma malha correspondente e assim permitir o cálculo de forças. A função genérica [`buildingslice`](@ref) tem métodos para fatiar o edifício. Neste nosso exemplo, vamos assumir que cada pavimento tem 3 m de altura.


```@example 3
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


### Calculando forças

Com o edifício discretizado, as forças podem facilmente ser calculadas com os objetos [`NodeInfo`](@ref) obtidos com [`buildsurface`](@ref) e [`buildingslice`](@ref) acima.

A partir de uma distribuição de pressão, a força agindo em cada face da superfície é dada por

$\vec{F} = -\int_S p \:d\vec{A}$

The moment is calculated by

$\vec{M} = -\int_S p \vec{r}\times d\vec{A}$

Using the discretization obtained above,

$\vec{F} = -\sum_{i=1}^{N_{taps}} p_i \cdot \vec{A}_i$

$\vec{M} = -\sum_{i=1}^{N_{taps}} p_i \cdot \vec{r_i} \times \vec{A}_i$

Repare que dado um vetor com todas as medições de pressão, a operação acima pode ser representada por uma multiplicação de matriz:

$\left\{\begin{matrix}F_x\\F_y\\F_z\\M_x\\M_y\\_Mz\end{matrix}\right\} = \left[ F_{matrix} \right] \cdot \left\{\begin{matrix}p_1\\p_2\\p_3\\\vdots\\p_{N_{taps}}\end{matrix}\right\}$

A matriz $\left[F_{matrix}\right] é esparsa. O número de linhas corresponde ao número de triângulos (o nós de maneira geral) e o número de colunas corresponde ao número de tomadas de pressão. Para calcular esta matriz para uma superfície discretizada, use o método [`addforcecontrib!`](@ref) ou [`forcematrix`](@ref). O método `forcematrix` alloca memória para a matriz e chama `addforcecontrib!` para calcular a contribuição da face de uma superfície. Os métodos são definidos assim pois pode haver contribuições de ambos os lados e estas contribuições são calculadas de maneira independente.

O primeiro argumento é o número de colunas que corresponde ao número de tomadas de pressão. O segundo argumento é a malha (`BuildingSurface`) e o terceiro argumento especifica quais componentes de força devem ser calculados:

 1. $F_x$
 2. $F_y$
 3. $F_z$
 4. $M_x$
 5. $M_y$
 6. $M_z$


```@example 3
# Lembre-se: tempos 240 tomadas de pressão!
Fbase = forcematrix(240, msh.nodes, (1,2,3,4,5,6); sgn=1, side=1, point=Point(0,0,0))
println("Dimensões de `Fbase`: $(size(Fbase))")
```

O argumento chave `sgn` multiplica cada elememento da matriz. No caso de pressões internas, este argumento geralmente é -1. Mas pode ser utilizado para levar em conta fatores de escala, conversão de unidade ou até mesmo velocidade de referência do vento. O argumento chave `side` especifica qual lado da superfície está contribuindo para a força. Em geral, `forcematrix` seria chamado com `sgn=1` e `side=1` para ter calcular a contribuição da face externa. Em seguida, `addforcecontrib!` é chamado com `sgn=-1` e `side=2` para adicionar a contribuição da pressão interna. Como neste primeiro exemplo existem apenas tomadas de pressão externa, apenas a função `forcematrix` será chamada com `sgn=1` e `side=1`.

Os momentos de cada triângulo são calculados em relação ao ponto especificado pelo argumento chave `point`.

No exemplo acima, a matriz de força construída calcula o carregamente na fundação do edifício.

### Calculando as forças em cada pavimento

Normalmente o calculista deseja a distribuição de carga ao longo da estrutura. No caso de edifícios altos, os programas de elementos finitos geralmente precisam das forças em cada pavimento. Se cada região da estrutura tem uma malha diferente, os métodos [`forcematrix`](@ref) e [`addforcecontrib!`](@ref) podem ser calculados de maneira independente. Mas esta é uma operação tão comum que foram implementados métodos para que calculam a matriz de força para vetores com diferentes superfícies. O uso é igual com um vetor de `BuildingSurface` usado no lugar de um único objeto `BuildingSurface`. Ainda existe um argumento chave `interleaved` que especifica em que sequência as forças são calculadas:


 * `interleaved=false`: cada componente da carga é numerada em sequência. Por exemplo, se o argumento `forces=(1,2,6)` a ordem das forças (cada linha da matriz de força) é  $F_{x,1}$, $F_{x,2}$, $\ldots$, $F_{x,N}$, $F_{y,1}$, $F_{y,2}$, $\ldots$, $F_{y,N}$, $M_{z,1}$, $M_{z,1}$, $\ldots$, $M_{z,N}$ onde $N$ é o número de pavimentos (mais especificamente o número de malhas).
 * `interleaved=false`: As cargas são numeradas primeiro por pavimento,  (ou malha), $F_{x,1}$, $F_{y,1}$, $M_{z,1}$, $F_{x,2}$, $F_{y,2}$, $M_{z,2}$, $\ldots$, $F_{x,N}$, $F_{y,N}$, $M_{z,N}$.


```@example 3
sl_nodes = [sl.nodes for sl in slices]
Fslices = forcematrix(240, sl_nodes, (1,2,6); sgn=1, side=1, point=Point(0,0,0))
println("Dimensions of `Fslices`: $(size(Fslices))")
```

## Um exemplo mais complexo

Este edifício tem todos os detalhes da primeira figura. O cilindro tem uma subdivisão que divide o cilindro em duas metades. Uma metade é isolada e não tem tomadas de pressão internal Mas a outra tem tanto tomadas de pressão externa quanto interna.

Este edifício tem duas superfícies:

 1. A superfície cilíndrica que tem duas seções:
    * Uma metade com tomadas de pressão tanto externas quanto internas
    * Uma metade com tomadas de pressão externas apenas.
 2. A superfície plana que divide o cilindro em duas metades e possui apenas tomadas de pressão externa.

### Definindo a geometria

Começaremos reproduzindo a geometria do edifício simples acima.

```@example 4
using Meshes, MeshViz
import GLMakie
using BuildingGeometry

H = 150.0 # Altura
D = 30.0  # Diametro
R = D/2   # Raio

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
```

### Definindo as tomadas de pressão externas e internas

Para a superfície 1, as tomadas de pressão externas são iguais ao caso do edifício simples acima. Mas temos que adicionar as tomadas de pressão interna em metade da superfície.

```@example 4
nr = 10 # Numero de linhas
dz = H/nr  # Altura de cada linha
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

# Nós internos da face 1
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


# Nós externos da face 2

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


### Discretizando o edifício

Agora vamos discretizar as superfícies do edifício. Este caso é mais complexo que o edifício simples acima. Primeiro que temos duas superfícies. A superfície cilíndrica externa e a divisão interna. Mas existe ainda outra dificuldade. Na superfície cilíndrica não existem tomadas de pressão internas em metade da superfície. Assim, esta superfície deveria ser decomposta em duas regiões.

Ainda existem outros aspectos que devem ser levados em consideração. O projetista quer o carregamento por pavimento mas na região aberta com tomadas internas e externas assim é interessante poder diferenciar as regiões de uma superfície. Para tratar de casos como este, [`NodeInfo`](@ref) tem um campo `tag` que pode associar um `Int` a cada lado da superfície. Assim, o edifício será decomposto como descrito a seguir:

 1. `face1` com 240 tomadas de pressão externa, tag externa 1
    1. Seção com 18 tomadas de pressão interna, tag interna 2
    2. Seção sem tomadas de pressão interna, tag interna 0
 2. `face2` com 9 tomadas de pressão externa e nehuma pressão interna, tag externo 3, tag interno 0 (o mesmo que para o tag interno da `face1`, seção 2)

As tomadas de pressão estão numeradas sequencialmente:

 1. Tomadas externas na superfície cilíndrica (1-240), coordenadas `epts1`
 2. Tomadas externas no plano interno (241-249), coordenadas `epts2`
 3. Tomadas internas em metada da superfície cilíndrica  (250-267)

A superfície do cilindro é caracterizada pelos triângulos no vetor `face1`. A primeira metade deste vetor (1-24) é a região com tomadas internas. Os triângulos 25-48 formam a região sem tomadas de pressão interna.


```@example 4

# Vamos discretizar a superfície cilíndrica 
msh1 = buildsurface(face1, # Tranalhando na primeira superfície
       		    [(points=epts1, tri=1:48, id=1:240, tag=1), # Face externa
		     (points=ipts1, tri=1:24, id=250:267, tag=2), # Face int, com tom
		     (tri=25:48, id=-1)]) # Face int, sem tom

# Esta chamada usa the `tag=0` and `nointid=-1` para a face interna
msh2 = buildsurface(face2, [(points=epts2, tri=1:2, id=241:249, tag=3)],
       			      nointid=-1, tag=0)

# Isto seria uma alternativa.
#msh2 = buildsurface(face2, [(points=epts2, tri=1:2, id=241:249, tag=3),
#      			     (tri=1:2, id=-1, tag=0)])

# Juntando as malhas

msh = mergemeshes(msh1, msh2)
```


Vamos ver algumas das faces

```@example 4
itags = nodetag.(msh.nodes, 2) # Internal tags

# We will first visualize the faces with no
# internal pressure taps
idx1 = itags .== 0

# Now let's checkout the faces with internal pressure taps
idx2 = .! idx1

# Lets view them:

let
   fig = GLMakie.Figure()
   ax1 = GLMakie.Axis3(fig[1,1], aspect=:data, title="No internal taps")
   viz!(tri2mesh(msh.tri[idx1]))

   ax2 = GLMakie.Axis3(fig[1,2], aspect=:data, title="Internal taps")
   viz!(tri2mesh(msh.tri[idx2]))
   
   GLMakie.save("figures/building3.png", GLMakie.current_figure());
end
```

![Faces with and without internal pressure taps](figures/building3.png)


### Forças

O mesmo que para o exemplo anterior. Com tags podemos especificar onde queremos
calcular as forças. Neste exemplo devemos lembrar das faces internas.


```@example 4
# Lembre-se, temos 267 tomadas de pressão!
Fbase = forcematrix(267, msh.nodes, (1,2,3,4,5,6); sgn=1, side=1, point=Point(0,0,0))

# Adicionar a contribuição das tomadas de pressão internas - tag=2
ii = nodetag.(msh.nodes, 2) .== 2 # Selecionando os nós com tag interna 2
addforcecontrib!(Fbase, msh.nodes[ii], (1,2,3,4,5,6); sgn=-1, side=2, point=Point(0,0,0))

println("Dimensões de `Fbase`: $(size(Fbase))")
```



## Visualizando os resultados

O pacote [`Makie`](https://docs.makie.org/stable/) é um poderoso sistema de visualização. Neste momento, `BuildingGeometry` usa o pacote [`Meshes.jl`](https://github.com/JuliaGeometry/Meshes.jl) para coisas geométricas ao invés do pacote `GeometryBasics.jl` usado por `Makie`. Como se pode ver nos exemplos acima, foi usado o pacote  [`MeshViz`](https://github.com/JuliaGeometry/MeshViz.jl) para visualizar as malhas. Este pacote faz uma ponte entre  `Meshes.jl`  e `Makie`.

Como exemplo, vamos tentar plotar na malha a função 

$f(x,y,z) = z \cdot \left[(R + x)^2 + 2(R+y)^2\right]$

```@example 4

function fun(p,R)
   x,y,z = coordinates(p)
   return z * ((R+x)^2 + 2*(R+y)^2)
end

u = fun.(msh.points, R)
smsh = tri2mesh(msh.tri)

data = meshdata(smsh, etable=(u=u,))
viz(data)
GLMakie.save("figures/building4.png", GLMakie.current_figure());
```

![Data visualization](figures/building4.png)


### [`WriteVTK.jl`](https://github.com/jipolanco/WriteVTK.jl)


O pacote [`WriteVTK.jl`](https://github.com/jipolanco/WriteVTK.jl) permite exportar resultados no formato VTK. Com arquivos neste formato, os resultados podem ser visualizados com ferramentas como  [Paraview](https://www.paraview.org/) ou [VisIt](https://visit-dav.github.io/visit-website/index.html).

