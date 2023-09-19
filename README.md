# The Robust Randomized Quasi Monte Carlo method, applications to integrating singular functions

This repository contains the code of the paper [The Robust Randomized Quasi Monte Carlo method, applications
to integrating singular functions](https://hal.science/hal-03631879) by E. Gobet M. Lerasle and D. Métivier.

The paper contains multiple figures. In all, the code $M$ is the number of realizations used to show the estimators' density. Since we are interested in the tails of this distribution, we typically need $M$ very large, which makes the simulation quite heavy. Some simulations lightweight and can run in less than an hour, they are in a Julia Pluto notebook.
Other heavy simulations are given a script meant to run in parallel on a Slurm or SGE cluster.

For all the Randomized Quasi Monte Carlo simulations, we use $M$ Nested Uniform Scramble of the Sobol' sequence[^1].

[^1]: The "true" Sobol' sequence where the first point starts at 0. Note that the [Sobol.jl](https://github.com/JuliaMath/Sobol.jl)@v.1.4 package used in [QuasiMonteCarlo.jl](https://github.com/SciML/QuasiMonteCarlo.jl) implements a truncated version. Instead, one can use the implementation in the (depreacted) [RandomizedQuasiMonteCarlo.jl](https://github.com/dmetivie/RandomizedQuasiMonteCarlo.jl) package or in [QMCGenerators.jl](https://github.com/alegresor/QMCGenerators.jl).

## Lightweight simulations

These simulations take less than an hour and can then be put into a Pluto Notebook.
The best way to visualize the notebooks is to open the `html` version of the notebook in a browser.
From there, you can "run or edit" with Julia (if installed on your machine) or with the free[^3] Binder cloud service.
The user is free to change parameters to go beyond what is shown in the paper.

- Figure 1 - Exactly reproducible in the self-contained Pluto notebook [`Figure_1.jl`](https://github.com/dmetivie/Robust-Randomized-Quasi-Monte-Carlo-paper-code/blob/5e702cf6dede8b50c8f6bb1c6faa961b04007044/short_simulation_notebooks/Figure_1.jl) or open the friendly [`html` version of `Figure_1.jl`](https://raw.githack.com/dmetivie/Robust-Randomized-Quasi-Monte-Carlo-paper-code/b8cd89c97883b4ebb1075fd1e605d0ec4a4af575/short_simulation_notebooks/Figure_1.html) in a web browser to see directly the resulting notebook. It runs for approximately 10 minutes with the paper settings of $M = 10^7$. This notebook only showcases robust estimators without any Quasi Monte Carlo.
  - Figure 1.a
  - Figure 1.b
- Figure 5 - Exactly reproducible in the Pluto notebook. The two following ``versions`` are exactly the same, the first one is just rendered as an `html` page.
  - The friendly notebook [version of `Figure_5.jl`](https://raw.githack.com/dmetivie/Robust-Randomized-Quasi-Monte-Carlo-paper-code/5e702cf6dede8b50c8f6bb1c6faa961b04007044/short_simulation_notebooks/Figure_5.html).
  - The script version of the notebook [`Figure_5.jl`](https://github.com/dmetivie/Robust-Randomized-Quasi-Monte-Carlo-paper-code/blob/5e702cf6dede8b50c8f6bb1c6faa961b04007044/short_simulation_notebooks/Figure_5.jl).
It runs for approximately 40 minutes with the paper settings of $M = 10^4$. This notebook shows the whole workflow to do multiple Robust Randomized Quasi Monte Carlo. Thanks to Pluto's capabilities, we only show the interesting and interactive pieces and hide[^2] the background functions.
- Figure 5 - Exactly reproducible in the Pluto notebook. The two following ``versions`` are exactly the same, the first one is just rendered as an `html` page.
  - The friendly [version of `Figure_5.jl`](https://raw.githack.com/dmetivie/Robust-Randomized-Quasi-Monte-Carlo-paper-code/5e702cf6dede8b50c8f6bb1c6faa961b04007044/short_simulation_notebooks/Figure_5.html).
  - The script version of the notebook [`Figure_5.jl`](https://github.com/dmetivie/Robust-Randomized-Quasi-Monte-Carlo-paper-code/blob/5e702cf6dede8b50c8f6bb1c6faa961b04007044/short_simulation_notebooks/Figure_5.jl)  
It runs for approximately 40 minutes with the paper settings of $M = 10^4$. This notebook shows the whole workflow to do multiple Robust Randomized Quasi Monte Carlo. Thanks to Pluto's capabilities, we only show the interesting and interactive pieces and hide[^2] the background functions.

[^2]: Hidden cells can be shown again by clicking on the eye symbol next to them. All the code is visible in the `.jl` script.
[^3]: This cloud service is (currently) free without registration. Hence, in one click, it loads the notebook and runs it. Setup might take some time.

## Heavy simulations

These simulations use high-dimensional RQMC simulations with $M = 10^4$. In order to run, one must typically use HPC/cluster computing. Unfortunately, Julia distributed computing is not very easily and efficiently reproducible, i.e., getting the same result regardless of the number of cores. It is not as simple as fixing a `seed` in serial computation.
Hence, here we provide the script that can be used to produce similar figures. Since $M\times n$ is quite high, the result should be very similar to the paper figure.

- Figure 2
  - Figure 2.a
  - Figure 2.b
- Figure 3
  - Figure 3.a
  - Figure 3.b
- Figure 4
- Figure 6

## Programming language and tooling used

### Julia

The [Julia programming language](https://julialang.org/) is used in this project.
One advantage of Julia is that it is both readable and fast (compiled).
There is no `C/C++/Fortan` wrapper.
All the code is Julia (and the packages used too). Hence, it is very easy to take a look at the randomization methods used in the paper, see [here](https://github.com/SciML/QuasiMonteCarlo.jl/blob/4cb665d27bed64871039b873e68060ef8cc374de/src/RandomizedQuasiMonteCarlo/scrambling_base_b.jl#L66) for example for Nested Uniform Scrambling.

Main packages used

- [QuasiMonteCarlo.jl](https://github.com/SciML/QuasiMonteCarlo.jl) developed by multiple authors. The documentation can be found [here](https://docs.sciml.ai/QuasiMonteCarlo/stable/). @dmetivie contributed to the Quasi Monte randomization method and docs of the package.
- [RobustMeans](https://github.com/dmetivie/RobustMeans.jl) developed by @dmetivie, one of the authors of the paper.

### Pluto notebooks

Most readers are probably familiar with Jupyter notebooks, Pluto's notebooks have plenty of noteworthy differences, see the [official website](https://plutojl.org/).
Here are some:

- Works only for Julia.
- Reactive, i.e., if you change one parameter somewhere in the notebook, all other parameters depending on function, plots, text, etc. will be automatically reevaluated to the new value. This allows one to quickly explore the effects of changing parameters, methods, etc.
- They are reproducible by nature, i.e., all the exact information on versioning used to generate the result is contained in the notebook.
- They are Julia files, i.e., one can open them as a simple and readable Julia script or as the notebook version.

### Clusters

Part of this work used the École Polytechnique [IDCS mesocentre](https://meso-ipp.gitlab.labos.polytechnique.fr/user_doc/) (Cholesky).
The final part of the work has been done using the [MESO@LR cluster](https://meso-lr.umontpellier.fr/documentation-utilisateurs/).
