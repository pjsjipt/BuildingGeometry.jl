
# Installation

The `BuildingGeometry` package is not registered so it should be installed directly:

From the REPL, in the package management mode,

```julia-repl
(tutorial) pkg> add https://github.com/pjsjipt/BuildingGeometry.jl
```

The package depends on ta few packages but for computing the Voronoi diagram it calls the Python package [scipy](https://https://scipy.org/). Probably a separate package should should be doing this but for now, your Python installation used by [`PyCall.jl`](https://github.com/JuliaPy/PyCall.jl) should be able to load scipy, in particular [`scipy.spatial.Voronoi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Voronoi.html)

In Linux, be sure to have the appropriate package installed. In Windows, once you add `BuildingGeometry`, it will add `PyCall` as a dependency and then you need to install the `scipy` package. The Julia package [`Conda.jl`](https://github.com/JuliaPy/Conda.jl) can help you with that.

If you want to use some other Python installed in your system, the [Anaconda distribution](https://www.anaconda.com/products/distribution) for instance, you will need to inform `PyCall.jl` which Python it should use. This can be done by setting the `PYTHON` environment variable to the path of the Python executable you want. After setting this env variable, don't forget to rebuild `PyCall`: in the package manager environment of the REPL, just type

```julia-repl
(tutorial) pkg> build PyCall
```

Now you should be ready to go...
