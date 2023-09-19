### A Pluto.jl notebook ###
# v0.19.27

#> custom_attrs = ["hide-enabled"]
#> 
#> [frontmatter]
#> Author = "David Métivier"
#> title = "Robust Quasi Monte Carlo Figure 5"
#> date = "2023-09-15"
#> tags = ["Quasi-MonteCarlo", "Robust-Statistic", "paper-with-code", "julia", "pluto-notebook"]
#> description = "This notebook shows the whole workflow to simulate and plot Figure 5 of the paper \"The Robust Randomized Quasi Monte Carlo method, applications to integrating singular functions\""

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 1cb37c25-2470-4094-8e09-d2ba44e8d68d
begin
    import Pkg
    # careful: this is _not_ a reproducible environment
    # activate the global environment
    Pkg.activate(".")
end


# ╔═╡ 45b56d47-11fb-43b8-9488-097244ef0919
begin
    using PlutoExtras
    using PlutoExtras.StructBondModule # This submodule has the features listed in this notebook
    using PlutoUI
    import PlutoUI: combine
    using HypertextLiteral
end

# ╔═╡ 2483fc8e-7fc9-4a45-89be-c2bb73fb10fe
using Printf

# ╔═╡ cfd51eb7-f484-4589-9af6-eccad8465b04
using OrderedCollections

# ╔═╡ 693fdf6b-884d-4502-b151-6ec3d4947741
using Random # for seeding

# ╔═╡ dbe828cf-79c6-4096-ab6f-db973c11ea59
using Distributions

# ╔═╡ a154d5f2-568e-4617-b972-e1a62b4c5081
using DelimitedFiles

# ╔═╡ adb19ce1-68e9-4840-8446-da36426e45b1
using Dates

# ╔═╡ d483b331-8dbb-4a25-9e3a-952681f065b8
using RobustMeans

# ╔═╡ 4f564974-2fbb-4151-bf29-b0f9b523ee45
using SpecialFunctions, HypergeometricFunctions

# ╔═╡ 7f75dd3f-b691-42ca-b0a4-d41f25e4fa48
using StatsPlots, LaTeXStrings

# ╔═╡ dc5b3731-1ebf-4cbb-a5f0-ceccdecc8ae9
using QuasiMonteCarlo, QMCGenerators

# ╔═╡ e79d3b3c-981f-4c3f-9fea-23fe544fb9f2
md"""
In this notebook, we focus on the scientific code and hide by default the code used to create the interactive features. This hidden code can be reveal in the notebook or seen directly in the script version.
"""

# ╔═╡ 0d7dc098-57b8-4960-bef1-f0d7fbf9bab3
md"""
# Packages
"""

# ╔═╡ 36e3cb04-bf9a-495f-ba2b-671d95a0ea94
# ╠═╡ custom_attrs = ["toc-hidden"]
md"""
## Pluto packages
"""

# ╔═╡ 4ce6403c-ae28-405e-af5f-3058f1bb4421
# html"""
# <style>
# 	@media screen {
# 		main {
# 			margin: 0 auto;
# 			max-width: 1500px;
#     		padding-left: max(283px, 10%);
#     		padding-right: max(383px, 10%); 
#             # 383px to accomodate TableOfContents(aside=true)
# 		}
# 	}
# </style>
# """

# ╔═╡ 7c29c5ec-a05c-450d-b108-d87abff57635
# ╠═╡ custom_attrs = ["toc-hidden"]
md"""
## Pluto nice looking stuff
"""

# ╔═╡ 1a3635dd-d283-4dfd-aa46-138a36166970
ExtendedTableOfContents(title="Robust Quasi Monte Carlo")

# ╔═╡ b0e1fbc2-eb39-425f-80f3-13f1ce64259d
struct Bypass
    content
end

# ╔═╡ 5441691c-077c-48a9-969c-024ed51da2a7
macro NTBondHack(desc, block)
    desc = esc(desc)
    :(@NTBond($Bypass($desc), $block))
end

# ╔═╡ f728941b-4280-4195-ab78-58bc8d39182a
Base.show(io::IO, f::Float64) = @printf(io, "%.2E", f) # scientific notation with two digits

# ╔═╡ 1391a0eb-8625-4a03-94ce-e6de6eef4c32
function Base.show(io::IO, ::MIME"text/javascript", o::Bypass)
    write(io, "`")
    Base.show(io, MIME"text/html"(), o.content)
    write(io, "`")
end

# ╔═╡ 2883121d-f293-4de3-a61f-a542caf59fe5
# using ProgressMeter

# ╔═╡ 3cdb1c90-569e-4d57-8b48-a0701966fb80
# ╠═╡ custom_attrs = ["toc-hidden"]
md"""
## Statistics packages
"""

# ╔═╡ 6928c3e7-6f14-47a3-843a-e1efdca23b17
Random.seed!(1234)

# ╔═╡ 5b78eba6-556a-4339-ad89-f00f9166939a
# ╠═╡ custom_attrs = ["toc-hidden"]
md"""
## Helpers
"""

# ╔═╡ 4c8bd31f-caac-494b-afaa-b947a877467c
function RobustMeans.mean(A::AbstractArray, δ::Real, Estimator::Catoni{<:Nothing}, kwargs...)
    n = length(A)
    α = RobustMeans.α_Catoni(δ, n, std(A))
    z = RobustMeans.Z_Estimator(α, RobustMeans.ψ_Catoni)
    return mean(A, z; kwargs...)
end

# ╔═╡ 0bc87c1d-70bf-4414-93be-5f3f25b88f95
function resample_inf!(QCM_dict, iterators, op)
    for (i, (key, val)) in enumerate(QCM_dict)
        a = 1
        while a > 0
            a = count(isinf, val)
            if a > 0
                iterators[i].count = a
                @info "$(string(key)) has $(a) inf"
                idx_inf = findall(isinf, val)
                val[idx_inf] = [op(X) for X in iterators[i]]
            end
        end
    end
end

# ╔═╡ f3cd9aef-b322-4eb4-9c48-fa59fd38e453
"""
    resample_inf(Upoints, iterators, op)
Resample inf points. Inf point can happen when evaluating singular function. In theory singular point are sampled with probability zero.
In practice, when working with finite bit expansion, it happens with nonzero probability.
Hence with just resample these point.
"""
function resample_inf(QCM_dict, iterators, op)
    QCM_dict_copy = copy(QCM_dict)
    resample_inf!(QCM_dict_copy, iterators, op)
    return QCM_dict_copy
end

# ╔═╡ d170d962-e228-489d-985e-73e040e59e1a
function RobustMeans.mean(A::AbstractArray, δ::Real, Estimator::Huber{<:Nothing}, kwargs...)
    n = length(A)
    α = RobustMeans.α_Huber(δ, n, std(A))
    z = RobustMeans.Z_Estimator(α, RobustMeans.ψ_Huber)
    return mean(A, z; kwargs...)
end

# ╔═╡ f338b66c-91d4-4858-b072-b1f1866fb781
"""
    computing_means(X::AbstractMatrix, estimators, δ)
Given X and an mean estimator compute the estimator for each column.
"""
function computing_means(X::AbstractMatrix, estimators, δ)
    results = OrderedDict{Symbol,Vector}()
    for (name, μ̂) in estimators
        results[name] = [mean(c, δ, μ̂) for c in eachcol(X)]
    end
    return results
end

# ╔═╡ 8391ed4c-770b-44a9-b740-1ac56b477744
function m1(a, b)
    return (1 / gamma(1 + a + b)) * 2^(-a - b) * gamma(1 / 2 + a) * (gamma(1 / 2 + b) +
                                                                     2^(1 / 2 + b) * gamma(1 + a + b) *
                                                                     pFq([1, 1 + a + b], [3 / 2 + a], -1) / gamma(3 / 2 + a))
end

# ╔═╡ 50e8e886-c872-453b-84e8-575849e380f4
function m2(a, b)
    return (2^(1 - 2 * a - 2 * b) * π * (csc(2 * b * π) + csc(2 * (a + b) * π)) * gamma(2 * a)) /
           (gamma(1 - 2 * b) * gamma(2 * (a + b))) +
           pFq([1 - 2 * a, 1 - 2 * a - 2 * b], [-2 * (-1 + a + b)], 1 / 2) / (-1 + 2 * a + 2 * b)
end

# ╔═╡ 630be454-326d-4c28-b578-6e0077b1d849
function g₂(x, a, b, m1_val, m2_val)
    denominator = sqrt(m2_val - m1_val^2)
    numerator = x^(b - 1 / 2) * abs(x - 1 / 2)^(a - 1 / 2) - m1_val
    return numerator / denominator
end

# ╔═╡ 4df97638-492a-49f5-835b-3af26ac2e181
function g₂(x, a, b)
    m1_val = m1(a, b)
    m2_val = m2(a, b)

    return g(x, a, b, m1_val, m2_val)
end

# ╔═╡ aec0ee12-a498-4696-83f1-0a55517aad0e
function g₂cos(x, a, b, c::Integer, m1_val, m2_val)
    denominator = sqrt(m2_val - m1_val^2)
    numerator = x^(b - 1 / 2) * abs(x - 1 / 2)^(a - 1 / 2) - m1_val

    return return (numerator + abs(cos(2π * (x - 0.1) * c)) - 2 / π) / denominator / 1.003660134355804
end

# ╔═╡ 8c3d40dc-0140-4a89-8fd6-d4828ddf2953
function short_name(name::AbstractString)
    s = split(name)
    return length(s) > 1 ? Symbol(first.(s)...) : Symbol(uppercase(string(first.(s, 2)...)))
end

# ╔═╡ 36abf44a-e974-4b35-887f-039a666dc50f
Base.@kwdef struct DigitalSeqB2G <: QuasiMonteCarlo.DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

# ╔═╡ 65af6cb0-5647-405b-8614-519cc5546f55
function QuasiMonteCarlo.sample(n, d, ::DigitalSeqB2G, T=Float64)
    return permutedims(Next(QMCGenerators.DigitalSeqB2G(d), n))
end

# ╔═╡ da42c3ed-f705-465a-ac43-d856d7789ad5
# ╠═╡ custom_attrs = ["toc-hidden"]
md"""
## Plotting
"""

# ╔═╡ 782c1627-f33a-48db-aa1b-07c16336798c
default(fontfamily="Computer Modern", linewidth=2, label=nothing)

# ╔═╡ a7dc3321-d33d-4b73-800d-4e6a2468c690
md"""
I just modify the whisker in boxplots so they are ajustable.
"""

# ╔═╡ b2fdc4ff-f785-44e9-93ac-02ce51f2e4f1
import StatsPlots: notch_width, _cycle

# ╔═╡ 35a539e0-23da-41ae-bb3a-29f1c0da6bbf
@recipe function f(
    ::Type{Val{:boxplot}},
    x,
    y,
    z;
    notch=false,
    whisker_range="0.95"::String,
    outliers=true,
    whisker_width=:half
)
    # if only y is provided, then x will be UnitRange 1:size(y,2)
    if typeof(x) <: AbstractRange
        if step(x) == first(x) == 1
            x = plotattributes[:series_plotindex]
        else
            x = [getindex(x, plotattributes[:series_plotindex])]
        end
    end
    xsegs, ysegs = Segments(), Segments()
    texts = String[]
    glabels = sort(collect(unique(x)))
    warning = false
    outliers_x, outliers_y = zeros(0), zeros(0)
    bw = plotattributes[:bar_width]
    isnothing(bw) && (bw = 0.8)
    @assert whisker_width == :match || whisker_width == :half || whisker_width >= 0 "whisker_width must be :match, :half, or a positive number"
    ww = whisker_width == :match ? bw :
         whisker_width == :half ? bw / 2 :
         whisker_width
    for (i, glabel) in enumerate(glabels)
        # filter y
        values = y[filter(i -> _cycle(x, i) == glabel, 1:length(y))]

        # compute quantiles
        q1, q2, q3, q4, q5 = quantile(values, range(0, stop=1, length=5))

        # notch
        n = notch_width(q2, q4, length(values))

        # warn on inverted notches?
        if notch && !warning && ((q2 > (q3 - n)) || (q4 < (q3 + n)))
            @warn("Boxplot's notch went outside hinges. Set notch to false.")
            warning = true # Show the warning only one time
        end

        # make the shape
        center = Plots.discrete_value!(plotattributes[:subplot][:xaxis], glabel)[1]
        hw = 0.5_cycle(bw, i) # Box width
        HW = 0.5_cycle(ww, i) # Whisker width
        l, m, r = center - hw, center, center + hw
        lw, rw = center - HW, center + HW

        # internal nodes for notches
        L, R = center - 0.5 * hw, center + 0.5 * hw

        # outliers
        qα = parse(Float64, whisker_range)
        @assert 0.5 < qα ≤ 1 "qα whisker_range must be in betweeen 0.5 and 1"
        if qα != 0.0  # if the range is 0.0, the whiskers will extend to the data
            limit_up = quantile(values, qα)
            limit_down = quantile(values, 1 - qα)

            inside = Float64[]
            for value in values
                if (value < limit_down) || (value > limit_up)
                    if outliers
                        push!(outliers_y, value)
                        push!(outliers_x, center)
                    end
                else
                    push!(inside, value)
                end
            end
            # change q1 and q5 to show outliers
            # using maximum and minimum values inside the limits
            q1, q5 = Plots.ignorenan_extrema(inside)
            q1, q5 = (min(q1, q2), max(q4, q5)) # whiskers cannot be inside the box
        end
        # Box
        push!(xsegs, m, lw, rw, m, m)       # lower T
        push!(ysegs, q1, q1, q1, q1, q2)    # lower T
        push!(
            texts,
            "Lower fence: $q1",
            "Lower fence: $q1",
            "Lower fence: $q1",
            "Lower fence: $q1",
            "Q1: $q2",
            "",
        )

        if notch
            push!(xsegs, r, r, R, L, l, l, r, r) # lower box
            push!(xsegs, r, r, l, l, L, R, r, r) # upper box

            push!(ysegs, q2, q3 - n, q3, q3, q3 - n, q2, q2, q3 - n) # lower box
            push!(
                texts,
                "Q1: $q2",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Q1: $q2",
                "Q1: $q2",
                "Median: $q3 ± $n",
                "",
            )

            push!(ysegs, q3 + n, q4, q4, q3 + n, q3, q3, q3 + n, q4) # upper box
            push!(
                texts,
                "Median: $q3 ± $n",
                "Q3: $q4",
                "Q3: $q4",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Q3: $q4",
                "",
            )
        else
            push!(xsegs, r, r, l, l, r, r)         # lower box
            push!(xsegs, r, l, l, r, r, m)         # upper box
            push!(ysegs, q2, q3, q3, q2, q2, q3)   # lower box
            push!(
                texts,
                "Q1: $q2",
                "Median: $q3",
                "Median: $q3",
                "Q1: $q2",
                "Q1: $q2",
                "Median: $q3",
                "",
            )
            push!(ysegs, q4, q4, q3, q3, q4, q4)   # upper box
            push!(texts, "Q3: $q4", "Q3: $q4", "Median: $q3", "Median: $q3", "Q3: $q4", "Q3: $q4", "")
        end

        push!(xsegs, m, lw, rw, m, m)             # upper T
        push!(ysegs, q5, q5, q5, q5, q4)          # upper T
        push!(
            texts,
            "Upper fence: $q5",
            "Upper fence: $q5",
            "Upper fence: $q5",
            "Upper fence: $q5",
            "Q3: $q4",
            "",
        )

    end

    if !Plots.isvertical(plotattributes)
        # We should draw the plot horizontally!
        xsegs, ysegs = ysegs, xsegs
        outliers_x, outliers_y = outliers_y, outliers_x

        # Now reset the orientation, so that the axes limits are set correctly.
        orientation := default(:orientation)
    end

    @series begin
        # To prevent linecolor equal to fillcolor (It makes the median visible)
        if plotattributes[:linecolor] == plotattributes[:fillcolor]
            plotattributes[:linecolor] = plotattributes[:markerstrokecolor]
        end
        primary := true
        seriestype := :shape
        x := xsegs.pts
        y := ysegs.pts
        ()
    end

    # Outliers
    if outliers && !isempty(outliers)
        @series begin
            primary := false
            seriestype := :scatter
            if get!(plotattributes, :markershape, :circle) == :none
                plotattributes[:markershape] = :circle
            end

            fillrange := nothing
            x := outliers_x
            y := outliers_y
            ()
        end
    end

    # Hover
    primary := false
    seriestype := :path
    marker := false
    if Plots.is_attr_supported(Plots.backend(), :hover)
        hover := texts
    end
    linewidth := 0
    x := xsegs.pts
    y := ysegs.pts
    ()
end

# ╔═╡ bb97a3bf-a917-4d33-9e8d-ed12bfeb368b
md"""
# Settings and Functions
"""

# ╔═╡ 4e58d175-5c38-4ca8-b10d-f16528a6ee7a
# ╠═╡ custom_attrs = ["toc-hidden"]
md"""
## Table bonds
"""

# ╔═╡ e4769581-bd49-43f2-955d-007ef8c36d07
bond_sample = @bind param_sample @NTBond "Sampling" begin
    M = (md"Number of realizations $M$", NumberField(1:1:10^9, default=10^4))
    d = (md"Dimension $d$", NumberField(1:20, default=7))
    base = (md"base", NumberField(2:1:17, default=2))
    m = (md"Exponent $m$", NumberField(1:20, default=11))
    n = (md"Sample $n$", NumberField(1:2^11, default=56))
    s = NumberField(0:30, default=0)
end

# ╔═╡ 296890ac-9282-498f-923a-6161935bb353
bond_func = @bind param_func @NTBond "Function" begin
    eps = (md"Singularity $\epsilon$", NumberField(0:1e-5:2, default=0.04))
    # c = (md"$\cos(2\pi*(x-0.1)c)$", NumberField(0:1:10^5, default = 20))
end

# ╔═╡ caf5d5fd-a668-45c3-8546-283439e7a329
bond_qmc = @bind param_qmc @NTBond "QMC" begin
    seq = (md" ", Select([DigitalSeqB2G => "Sobol_Cs", SobolSample => "Sobol wo 0", FaureSample => "Faure", HaltonSample => "Halton"]))
end

# ╔═╡ 7b715692-6990-4751-a29f-ac0fdd82de38
bond1 = BondTable([
        bond_sample,
        bond_func,
        bond_qmc
]; description="Parameters") |> show_output_when_hidden

# ╔═╡ d6365e9f-b1fe-4194-a7f0-608ba9392296
md"""
## Functions
"""

# ╔═╡ c83c3717-7edd-4f29-8cba-67c48b7aa3d5
md"""
Syntetic $f$ from [Chapter 15.9](https://artowen.su.domains/mc/practicalqmc.pdf) from Art B. Owen
"""

# ╔═╡ 0440acca-0731-4db3-9f50-017aa95c461e
md"""
```math
F(x) = F_{\beta, G}(x) = \prod_{k=1}^d (1+\beta_k G_k(x_k))
```

```math
\int_0^1 G_k(x) \mathrm{d}x = 0, \quad \int_0^1 G_k^2(x) \mathrm{d}x = 1
```

```math
\mu = 1, \quad \sigma^2 = \prod_{k=1}^d (1+\beta_k^2)-1
```
"""

# ╔═╡ ea7f9740-6795-4744-9b83-25466615ad0d
F(x::AbstractArray, G, β) = prod(1 + β[i] * G(xᵢ) for (i, xᵢ) in enumerate(x))

# ╔═╡ 4328c5fc-9fe9-48a0-a178-1706b4bb7792
σ_F(β) = sqrt(prod(1 + βₖ^2 for βₖ in β) - 1)

# ╔═╡ 2bdea7eb-8c96-45c4-adfd-cf090acaa0b8
md"""
We choose our function and parameters such that $\mu = 1$ and $\sigma \approx 1/2$ for all dimension $d$ (so that comparison is easier).
"""

# ╔═╡ 8c700188-3c4e-4107-aeda-8c5a2c73d898
μ_exact = 1

# ╔═╡ de21e3e4-8a36-4f51-8303-b11987293b98
md"""
## Parameters
"""

# ╔═╡ 268ed34a-f6f1-4624-bd63-d3acfb668dbc
m = param_sample.m - param_sample.s

# ╔═╡ 9143887c-23df-4d04-a43b-cdd1b92f585d
base = param_sample.base

# ╔═╡ e8908310-4fdb-420e-ba44-5278de904b53
N = base^m

# ╔═╡ ee0207ac-f3bd-4be8-888a-9d2e4bc116c9
pad = 32

# ╔═╡ 1e3fc402-802d-410c-b669-61ee567a75f0
n = param_sample.n * (base^param_sample.s)

# ╔═╡ b5b8217d-1902-410a-ad14-4501a899c9fe
M = param_sample.M

# ╔═╡ 1075c82a-d6db-45c9-b2b0-3564ee5dda67
const ϵ = param_func.eps

# ╔═╡ ebb3f903-4c62-4a11-b43b-18eca8b20184
md"""
The choice bellow ensure a fixed variance through when changing dimension (for fair comparison)
"""

# ╔═╡ d6d41089-39f1-477a-a5f5-df399eb3ba2d
β_values = [0.346574, 0.286243, 0.2615, 0.24686, 0.236718, 0.229051, 0.222928, 0.217849, 0.213519, 0.209752, 0.206421, 0.203437, 0.200736, 0.198269, 0.196, 0.1939, 0.191946, 0.190118, 0.188402, 0.186784]

# ╔═╡ 86d9be9c-408c-4e32-a128-bcb38af1c8bd
md"""
The dimension of integration $d=$
"""

# ╔═╡ 8c9365f6-c155-47de-94a5-007f2356bbb1
d = param_sample.d

# ╔═╡ 45d939a2-e47b-4660-bcbc-9eaf87b41de2
β = β_values[d] # could be codded as a function?

# ╔═╡ 126b0be0-10e8-46fb-b872-13eb3c5155f4
# β = 0.217849; d = 8
# # β = 0.209752; d = 10
# # β = 0.200736;    d = 13
# # β = 0.198269; d = 14
# # β = 0.196; d = 15
# # β = 0.1939; d = 16
# # β = 0.191946; d = 17
# # β = 0.179855 # d = 25
# # β = 0.169671 # d = 35

# ╔═╡ 9b9888c4-9b98-4046-8bf7-80a699b6da40
md"""
## Operations
"""

# ╔═╡ ec7239b1-79a8-454e-af57-f95334a73670
const βₛ = β ./ [log(1 + s) for s in 1:d]

# ╔═╡ 95f14109-bce7-469d-8f23-e64f0105787f
σ_exact = σ_F(βₛ)

# ╔═╡ dd50a850-20b3-49d9-804e-c8ea54d868a1
const a1 = m1(ϵ, ϵ)

# ╔═╡ 7f27cd1e-3186-4fc1-9f57-77680ed74a16
const a2 = m2(ϵ, ϵ)

# ╔═╡ bb2fb72e-fdad-40b5-819f-d81e6cd4aaca
const c = 20

# ╔═╡ 9f9c82d6-8965-4f44-b2af-efef5dfe5873
md"""
## Integrand
"""

# ╔═╡ a0e861e4-55f8-47c1-8d55-585f36b5b298
"""
	G_2sing(x)
Function with singularities at `x = (0,0,..., 0)` and `x = (1/2,1/2,..., 1/2)`.
"""
G_2sing(x) = g₂(x, ϵ, ϵ, a1, a2)

# ╔═╡ 78436d66-848a-453c-9917-b3a432e289e7
"""
	G_2sing_cos(x)
Function with singularities at `x = (0,0,..., 0)` and `x = (1/2,1/2,..., 1/2)` + a cos term.
"""
G_2sing_cos(x) = g₂cos(x, ϵ, ϵ, c, a1, a2)

# ╔═╡ be61b5ce-5909-4622-ad50-f81304c50196
"""
	G_1sing(x)
Function with a singularity at `x = (0,0,..., 0)`.
"""
G_1sing(x) = (sqrt(2ϵ) * (1 + 2ϵ)) / abs(1 - 2ϵ) * (x^(-1 / 2 + ϵ) - 2 / (1 + 2ϵ))

# ╔═╡ e51eba94-146b-4015-87d4-3284c998e7c4
md"""
Choose you integrand. In the paper `G_2sing` is used with the default parameters of this notebook to produce Figure 6.
"""

# ╔═╡ 7de68a6c-310a-4e48-af4f-39a096990c9f
@bind G Select([G_1sing, G_2sing, G_2sing_cos], default=G_2sing)

# ╔═╡ aae20b8a-8c92-4d1d-b89c-e8d6ecb4ac85
f(x) = F(x, G, βₛ)

# ╔═╡ 1ae2de43-59d7-4944-9dc7-e4f7c6a1cf2a
operation(f, X) = mean(f(c) for c in eachcol(X))

# ╔═╡ 242f166a-eb3e-47ac-a191-5b990c152f40
md"""
# Quasi Monte Carlo
"""

# ╔═╡ 744c350b-da55-4aa1-9c4d-6cbee45b324f
md"""
## QMC Sequence to be randomized
"""

# ╔═╡ 8a6e0b0c-d93c-4b13-920d-a68154b54a69
md"""
The default choice `DigitalSeqB2G` produce the correct Sobol' sequence, while the other Sobol truncate the first $(0,...,0)$ point and move other points.
"""

# ╔═╡ 422eb07b-9a17-4b63-83fd-ef7ff3b7187d
seq = param_qmc.seq

# ╔═╡ 02e689c0-3f1a-445a-8d11-f45b92d85c9e
rand_names = ["Monte Carlo", "Nestesd Uniform Scrambling", "Linear Matrix Scrambling", "Digital Shift", "Shift"]

# ╔═╡ bdfeb5c4-3792-4b0a-a87c-27213780833c
md"""
## Randomization
"""

# ╔═╡ f7b45bc6-4f41-4fc9-9ec7-0c35f032f03c
@bind rqmc_idx MultiCheckBox([i => n for (i, n) in enumerate(rand_names)], select_all=true, orientation=:column, default=[i for i in 1:2])

# ╔═╡ 4530f2c5-e3ee-48c1-9bc1-0676b7ece06a
md"""
## Loop
"""

# ╔═╡ 792b4553-9896-46f3-9c3c-71e077a77e9c
out_type = Float32

# ╔═╡ 1bf4249f-cfb7-43a9-af9a-360d86409b38
meth(b, pad, seq) = [RandomSample(); [seq(R=RQMC_meth(base=b, pad=pad)) for RQMC_meth in [OwenScramble, MatousekScramble, DigitalShift]]; Shift()]

# ╔═╡ 0a4fb17f-6514-429c-a861-7fd01e3dfc48
rqmc_methods = meth(Int32(base), Int32(pad), seq)[rqmc_idx]

# ╔═╡ b41e2083-5a66-43c7-84a0-363cc242615f
iterators = [DesignMatrix(N, d, algorithm, M * n) for algorithm in rqmc_methods]

# ╔═╡ a0651193-f621-4ee7-9777-f3209140db01
function run_long_randomization(iterators, rqmc_idx, rand_names)
    μ_QMC = OrderedDict{Symbol,Vector{out_type}}()

    for (i, iterator) in enumerate(iterators)
        name = rand_names[rqmc_idx[i]]
        println(name)
        sn = short_name(name)
        @time μ_QMC[sn] =
            map(iterator) do X
                mean(f(c) for c in eachcol(X))
            end
    end
    return μ_QMC
end

# ╔═╡ 3e7d3c18-c299-4bf6-a160-d518e01fba20
md"""
!!! warning "☡ Big computation happening here!"
	For `M = 100` this takes about 21.9 s, hence for `M = 10_000` as in the paper, expect around 2190s≈37min.
"""

# ╔═╡ 519ad352-0774-48ed-a5c9-702122a8d2ee
μ_QMC_unsafe = run_long_randomization(iterators, rqmc_idx, rand_names)

# ╔═╡ 6d17e797-5528-4cf2-966e-78382744dac6
md"""
## Resample `Inf`
"""

# ╔═╡ 2923f6be-73f1-4196-a935-06054289a1f8
md"""
When 0<ϵ<1/2, `f` has a singularity. Hence it can happen that our evaluation is `Inf`. In theory this happens with probability 0.
In practice, QMC we use 32 bits encoding for Floats, hence with a non zero probability after a large number of randomization we might encouter some singularities.
"""

# ╔═╡ ef2c07d4-9246-401c-a444-647a748f42df
μ_QMC = resample_inf(μ_QMC_unsafe, iterators, x -> operation(f, x))

# ╔═╡ 2bd3fe05-a7d3-4a87-821a-7c95a76e890f
md"""
# Robust Estimators
"""

# ╔═╡ 0dce832d-b462-4bc7-afb4-4643c966234e
md"""
## Selecting estimator, uncertainty $\delta$
"""

# ╔═╡ af759074-66fb-4bf2-b272-d8450c0848c7
md"""
Uncertainty level $\delta$.
This choice of δ ensures that there are no problems when grouping data into blocks
"""

# ╔═╡ 0c3f8914-402a-4add-8edf-8e9621af823b
δ = 3exp(-8) # ≈ 0.005

# ╔═╡ eb369c80-9456-4be1-9792-3fc4b3debf4f
robust_names = ["Empirical Mean", "Huber", "Median of Mean", "Catoni", "Lee Valiant", "Minsker Ndaoud", "Trimmed Mean"]

# ╔═╡ 0279854b-1851-4921-aa83-9ab86d88120a
md"""
Select the estimator you want to test.
"""

# ╔═╡ 32478a4f-7775-479d-a904-b56a3df3c5fc
md"""
!!! note
	Trimmed Mean require $n>(32/3)/\log(8/\delta)$, which for our default $\delta$ is $n > 96$. Hence, we do not use it here. But if you change to a larger $n$ you will be able to enable Trimmed Mean without errors.
"""

# ╔═╡ 77edddb4-4dff-4e1c-b435-f41d810de633
md"""
Select you Robust estiamtor
"""

# ╔═╡ 5e8ffc0d-a485-46cc-9319-f09169dc0ec1
@bind estimator_idx MultiCheckBox([i => n for (i, n) in enumerate(robust_names)], select_all=true, orientation=:column, default=[i for i in 1:(length(robust_names)-1)])

# ╔═╡ bb3b6494-145e-4acd-932e-fd50c929412d
md"""
Minsker Ndaoud hyperparameter
"""

# ╔═╡ e26d4c6e-da1f-43be-9cd8-cdca170424bd
pMN = 1

# ╔═╡ f11ea168-9b81-4658-8311-f3b3c5079b46
estim() = [EmpiricalMean(), Huber(nothing), MedianOfMean(0), Catoni(nothing), LeeValiant(), MinskerNdaoud(pMN), TrimmedMean(0, 0)]

# ╔═╡ 33346008-fde2-4ddc-ac8a-2be48b3ce923
md"""
## Computing the mean estimators $\hat{\mu}_n$
"""

# ╔═╡ abf5cf90-d2d1-4acb-815a-82e3094e6469
estimators = OrderedDict(short_name(robust_names[i]) => estim()[i] for i in estimator_idx)

# ╔═╡ ec726b24-91e0-4111-9e85-171ab6bb239b
md"""
!!! warning "Computation here"
"""

# ╔═╡ 60d9926c-df4d-4300-8dd8-82e6d384fa3b
results = OrderedDict{Symbol,OrderedDict}(name => computing_means(Float64.(reshape(U, n, M)), estimators, δ) for (name, U) in μ_QMC)

# ╔═╡ a0973218-bba8-4959-8892-98e0a1762bd4
md"""
# Plot
"""

# ╔═╡ 5c630246-91b5-4841-9fe3-312f13406602
md"""
## Box plot result
"""

# ╔═╡ fde89e74-822b-4892-a11e-de300d1c4d96
just_line = false

# ╔═╡ 352cb249-eb69-40d3-ad57-8b84fdce9850
li = just_line ? "line" : "box"

# ╔═╡ 67d02b48-9a44-40c4-ad91-ca4c92d48cc2
seqs_name = length(μ_QMC) == 2 ? ["MC", "RQMC"] : collect(keys(μ_QMC))

# ╔═╡ 9e61b2ca-5e28-42b3-b251-5648a2eb34e8
sequence = keys(μ_QMC)

# ╔═╡ 09d8c0d2-ea31-4720-b74f-537c0b8c2424
sequence

# ╔═╡ 59d66617-0b42-474b-881d-9244f62d3096
estimator = [:EM, :LV, setdiff(collect(keys(estimators)), [:EM, :LV])...]

# ╔═╡ c898f2d0-b57b-4827-9a0f-43fd0dc0f7c9
begin
    x_coord = zeros(Int, 0)
    y_coord = zeros(0)
    z_coord = zeros(Int, 0)

    plot(thickness_scaling=2)
    for (seq_num, s) in enumerate(sequence)
        for (kk, estimate) in enumerate(estimator)
            val = results[s][estimate] / μ_exact
            ΔM = abs.(val .- 1)
            println("$s $estimate 𝔼(μ̂) - μ = ", mean(val .- 1), " σ(μ̂) = ", std(val))
            println("RMSE(μ̂) = ", sqrt(var(val) + mean(val .- 1)^2))
            append!(y_coord, ΔM)
            append!(x_coord, fill(kk, length(ΔM)))
            append!(z_coord, fill(seq_num, length(ΔM)))
            if just_line == true
                if kk == 1
                    plot!([kk + (seq_num - 3 / 2) * 0.3 - 0.1, kk + (seq_num - 3 / 2) * 0.3 + 0.1], [quantile(ΔM, 1 - δ), quantile(ΔM, 1 - δ)], c=seq_num, label=seqs_name[seq_num])
                else
                    plot!([kk + (seq_num - 3 / 2) * 0.3 - 0.1, kk + (seq_num - 3 / 2) * 0.3 + 0.1], [quantile(ΔM, 1 - δ), quantile(ΔM, 1 - δ)], c=seq_num)
                end
                plot!([kk + (seq_num - 3 / 2) * 0.3, kk + (seq_num - 3 / 2) * 0.3], [quantile(ΔM, 0), quantile(ΔM, 1 - δ)], c=seq_num, label=:none, size=(1000, 500), alpha=1)
                plot!([kk + (seq_num - 3 / 2) * 0.3 - 0.1, kk + (seq_num - 3 / 2) * 0.3 + 0.1], [quantile(ΔM, 0), quantile(ΔM, 0)], c=seq_num)
                plot!([kk + (seq_num - 3 / 2) * 0.3 - 0.1, kk + (seq_num - 3 / 2) * 0.3 + 0.1], [quantile(ΔM, 1 - δ), quantile(ΔM, 1 - δ)], c=seq_num)

            end
        end
    end
    hline!([quantile(Normal(), 1 - δ / 2) / sqrt(n) * std(μ_QMC[:MC]) / μ_exact], c=1, label=:none, s=:dot)
    hline!([quantile(Normal(), 1 - δ / 2) / sqrt(n) * std(μ_QMC[:NUS]) / μ_exact], c=2, label=:none, s=:dash)
    if just_line == false
        groupedboxplot!(x_coord, y_coord, group=z_coord, label=permutedims(seqs_name), size=(1200, 600), alpha=1, c=permutedims(1:length(sequence)), ms=1.8,
            whisker_range="$(round(1-δ, digits = 6))"
        )
    end
    xticks!(1:length(estimator), string.(estimator), background_color_legend=nothing)
    yaxis!(:log10, yminorticks=9, minorgrid=:y, legend=(0.35, 0.93), minorgridlinewidth=1.0)
    ylims!((0.999991e-4, 1.01e-1))
    yticks!(10.0 .^ (-9:-0))
    ylabel!(L"|\hat{\mu}_{N,n}-\mu|/\mu", tickfonthalign=:center)
end

# ╔═╡ bbe9940c-1626-496b-80fa-00dfab30a031
md"""
We display the Root Mean Square error as well as the bias and variance. This shows that not only robust estimator are better than EM for the tails but are better at RMSE (not at bias only though).
"""

# ╔═╡ 20cd6866-1508-4f00-956e-c64337207eca
md"""
## Save figure
"""

# ╔═╡ d7831b42-a57b-4f34-ae0a-61f00e9f9c1c
begin # save the figure and crop it to avoid white spaces appearing in gr()
    save_name = "box_plot_$(G)_m_$(m)_n_$(n)_d_$(d)_esp_$(round(ϵ, digits = 2))_delta_$(round(δ, digits = 3))_$(string(estimator...))_M_$(M)_$(li)_$(today())"
    save_name = replace(save_name, "." => "p")
    savefig(string(save_name, ".pdf"))
    run(`pdfcrop $(string(save_name,".pdf"))`) # Petit délire pour croper proprement la figure 
    mv(string(save_name, "-crop", ".pdf"), string(save_name, ".pdf"), force=true)
end


# ╔═╡ 6a2bb211-e55e-4973-bff2-310085340e70
md"""
## Estimator distribution
"""

# ╔═╡ afd12f94-74fd-4f6f-a7b9-3aaf6a3e5ffe
md"""
Randomization method
"""

# ╔═╡ ca0483d3-fb9f-4035-a5dd-a6ab1da93605
@bind qmc_seq Select(collect(keys(μ_QMC)) .=> rand_names[rqmc_idx], default=collect(keys(μ_QMC))[2])

# ╔═╡ cc20226d-2d2e-4d5a-9936-d725e591eb45
begin
    fr = true # focus on central region ?
    framebox = fr ? "zoom" : "all"
    plot(thickness_scaling=2, size=(1000, 600))
    plot!(Normal(), label=L"\mathcal{N}(0,1)", c=:black, alpha=0.6)
    for (i, (name, s)) in enumerate(results[qmc_seq])
        W = √(n * N) * (s .- μ_exact) / σ_exact
        density!(W, alpha=0.6, norm=:pdf, label=string(name), c=i)
        # vline!([quantile(W, 1-δ)], s = :dot, c = i)
    end

    vline!([0], label=:none, c=:black, lw=1, alpha=0.9)
    yaxis!(:log10, yminorticks=9, minorgrid=:y, legend=:topright, minorgridlinewidth=1.2)
    ylims!((1 / M * 10, 2))
    xlabel!(L"\sqrt{n}(\hat{\mu}_n-\mu)/\sigma", tickfonthalign=:center)
    ylabel!("PDF")
    # if fr
    #     xlims!((-5, 10))
    # end
    # yticks!(10.0 .^ (-7:-0))
end

# ╔═╡ Cell order:
# ╟─0d7dc098-57b8-4960-bef1-f0d7fbf9bab3
# ╟─e79d3b3c-981f-4c3f-9fea-23fe544fb9f2
# ╟─36e3cb04-bf9a-495f-ba2b-671d95a0ea94
# ╠═1cb37c25-2470-4094-8e09-d2ba44e8d68d
# ╠═4ce6403c-ae28-405e-af5f-3058f1bb4421
# ╟─7c29c5ec-a05c-450d-b108-d87abff57635
# ╠═1a3635dd-d283-4dfd-aa46-138a36166970
# ╠═45b56d47-11fb-43b8-9488-097244ef0919
# ╟─b0e1fbc2-eb39-425f-80f3-13f1ce64259d
# ╟─5441691c-077c-48a9-969c-024ed51da2a7
# ╟─1391a0eb-8625-4a03-94ce-e6de6eef4c32
# ╠═2483fc8e-7fc9-4a45-89be-c2bb73fb10fe
# ╠═f728941b-4280-4195-ab78-58bc8d39182a
# ╠═2883121d-f293-4de3-a61f-a542caf59fe5
# ╟─3cdb1c90-569e-4d57-8b48-a0701966fb80
# ╠═cfd51eb7-f484-4589-9af6-eccad8465b04
# ╠═693fdf6b-884d-4502-b151-6ec3d4947741
# ╠═dbe828cf-79c6-4096-ab6f-db973c11ea59
# ╠═a154d5f2-568e-4617-b972-e1a62b4c5081
# ╠═adb19ce1-68e9-4840-8446-da36426e45b1
# ╠═d483b331-8dbb-4a25-9e3a-952681f065b8
# ╠═4f564974-2fbb-4151-bf29-b0f9b523ee45
# ╠═7f75dd3f-b691-42ca-b0a4-d41f25e4fa48
# ╠═dc5b3731-1ebf-4cbb-a5f0-ceccdecc8ae9
# ╠═6928c3e7-6f14-47a3-843a-e1efdca23b17
# ╟─5b78eba6-556a-4339-ad89-f00f9166939a
# ╟─f338b66c-91d4-4858-b072-b1f1866fb781
# ╟─4c8bd31f-caac-494b-afaa-b947a877467c
# ╟─0bc87c1d-70bf-4414-93be-5f3f25b88f95
# ╟─f3cd9aef-b322-4eb4-9c48-fa59fd38e453
# ╟─d170d962-e228-489d-985e-73e040e59e1a
# ╟─8391ed4c-770b-44a9-b740-1ac56b477744
# ╟─50e8e886-c872-453b-84e8-575849e380f4
# ╠═630be454-326d-4c28-b578-6e0077b1d849
# ╠═4df97638-492a-49f5-835b-3af26ac2e181
# ╠═aec0ee12-a498-4696-83f1-0a55517aad0e
# ╟─8c3d40dc-0140-4a89-8fd6-d4828ddf2953
# ╟─36abf44a-e974-4b35-887f-039a666dc50f
# ╟─65af6cb0-5647-405b-8614-519cc5546f55
# ╟─da42c3ed-f705-465a-ac43-d856d7789ad5
# ╠═782c1627-f33a-48db-aa1b-07c16336798c
# ╟─a7dc3321-d33d-4b73-800d-4e6a2468c690
# ╠═b2fdc4ff-f785-44e9-93ac-02ce51f2e4f1
# ╟─35a539e0-23da-41ae-bb3a-29f1c0da6bbf
# ╟─bb97a3bf-a917-4d33-9e8d-ed12bfeb368b
# ╟─4e58d175-5c38-4ca8-b10d-f16528a6ee7a
# ╠═7b715692-6990-4751-a29f-ac0fdd82de38
# ╟─e4769581-bd49-43f2-955d-007ef8c36d07
# ╟─296890ac-9282-498f-923a-6161935bb353
# ╟─caf5d5fd-a668-45c3-8546-283439e7a329
# ╟─d6365e9f-b1fe-4194-a7f0-608ba9392296
# ╟─c83c3717-7edd-4f29-8cba-67c48b7aa3d5
# ╟─0440acca-0731-4db3-9f50-017aa95c461e
# ╠═ea7f9740-6795-4744-9b83-25466615ad0d
# ╠═4328c5fc-9fe9-48a0-a178-1706b4bb7792
# ╟─2bdea7eb-8c96-45c4-adfd-cf090acaa0b8
# ╟─8c700188-3c4e-4107-aeda-8c5a2c73d898
# ╟─95f14109-bce7-469d-8f23-e64f0105787f
# ╟─de21e3e4-8a36-4f51-8303-b11987293b98
# ╟─268ed34a-f6f1-4624-bd63-d3acfb668dbc
# ╟─9143887c-23df-4d04-a43b-cdd1b92f585d
# ╟─e8908310-4fdb-420e-ba44-5278de904b53
# ╟─ee0207ac-f3bd-4be8-888a-9d2e4bc116c9
# ╟─1e3fc402-802d-410c-b669-61ee567a75f0
# ╟─b5b8217d-1902-410a-ad14-4501a899c9fe
# ╟─1075c82a-d6db-45c9-b2b0-3564ee5dda67
# ╟─ebb3f903-4c62-4a11-b43b-18eca8b20184
# ╟─d6d41089-39f1-477a-a5f5-df399eb3ba2d
# ╟─86d9be9c-408c-4e32-a128-bcb38af1c8bd
# ╟─8c9365f6-c155-47de-94a5-007f2356bbb1
# ╟─45d939a2-e47b-4660-bcbc-9eaf87b41de2
# ╟─126b0be0-10e8-46fb-b872-13eb3c5155f4
# ╟─9b9888c4-9b98-4046-8bf7-80a699b6da40
# ╟─ec7239b1-79a8-454e-af57-f95334a73670
# ╟─dd50a850-20b3-49d9-804e-c8ea54d868a1
# ╟─7f27cd1e-3186-4fc1-9f57-77680ed74a16
# ╟─bb2fb72e-fdad-40b5-819f-d81e6cd4aaca
# ╟─9f9c82d6-8965-4f44-b2af-efef5dfe5873
# ╟─a0e861e4-55f8-47c1-8d55-585f36b5b298
# ╟─78436d66-848a-453c-9917-b3a432e289e7
# ╟─be61b5ce-5909-4622-ad50-f81304c50196
# ╟─e51eba94-146b-4015-87d4-3284c998e7c4
# ╟─7de68a6c-310a-4e48-af4f-39a096990c9f
# ╟─aae20b8a-8c92-4d1d-b89c-e8d6ecb4ac85
# ╟─1ae2de43-59d7-4944-9dc7-e4f7c6a1cf2a
# ╟─242f166a-eb3e-47ac-a191-5b990c152f40
# ╟─744c350b-da55-4aa1-9c4d-6cbee45b324f
# ╟─8a6e0b0c-d93c-4b13-920d-a68154b54a69
# ╠═422eb07b-9a17-4b63-83fd-ef7ff3b7187d
# ╟─02e689c0-3f1a-445a-8d11-f45b92d85c9e
# ╟─bdfeb5c4-3792-4b0a-a87c-27213780833c
# ╟─f7b45bc6-4f41-4fc9-9ec7-0c35f032f03c
# ╟─4530f2c5-e3ee-48c1-9bc1-0676b7ece06a
# ╠═792b4553-9896-46f3-9c3c-71e077a77e9c
# ╠═1bf4249f-cfb7-43a9-af9a-360d86409b38
# ╟─0a4fb17f-6514-429c-a861-7fd01e3dfc48
# ╠═b41e2083-5a66-43c7-84a0-363cc242615f
# ╟─a0651193-f621-4ee7-9777-f3209140db01
# ╟─3e7d3c18-c299-4bf6-a160-d518e01fba20
# ╠═519ad352-0774-48ed-a5c9-702122a8d2ee
# ╟─6d17e797-5528-4cf2-966e-78382744dac6
# ╟─2923f6be-73f1-4196-a935-06054289a1f8
# ╠═ef2c07d4-9246-401c-a444-647a748f42df
# ╟─2bd3fe05-a7d3-4a87-821a-7c95a76e890f
# ╟─0dce832d-b462-4bc7-afb4-4643c966234e
# ╟─af759074-66fb-4bf2-b272-d8450c0848c7
# ╠═0c3f8914-402a-4add-8edf-8e9621af823b
# ╟─eb369c80-9456-4be1-9792-3fc4b3debf4f
# ╟─0279854b-1851-4921-aa83-9ab86d88120a
# ╟─32478a4f-7775-479d-a904-b56a3df3c5fc
# ╟─77edddb4-4dff-4e1c-b435-f41d810de633
# ╟─5e8ffc0d-a485-46cc-9319-f09169dc0ec1
# ╠═f11ea168-9b81-4658-8311-f3b3c5079b46
# ╟─bb3b6494-145e-4acd-932e-fd50c929412d
# ╟─e26d4c6e-da1f-43be-9cd8-cdca170424bd
# ╟─33346008-fde2-4ddc-ac8a-2be48b3ce923
# ╟─abf5cf90-d2d1-4acb-815a-82e3094e6469
# ╟─ec726b24-91e0-4111-9e85-171ab6bb239b
# ╠═60d9926c-df4d-4300-8dd8-82e6d384fa3b
# ╟─a0973218-bba8-4959-8892-98e0a1762bd4
# ╟─5c630246-91b5-4841-9fe3-312f13406602
# ╟─fde89e74-822b-4892-a11e-de300d1c4d96
# ╟─352cb249-eb69-40d3-ad57-8b84fdce9850
# ╟─67d02b48-9a44-40c4-ad91-ca4c92d48cc2
# ╟─9e61b2ca-5e28-42b3-b251-5648a2eb34e8
# ╟─09d8c0d2-ea31-4720-b74f-537c0b8c2424
# ╟─59d66617-0b42-474b-881d-9244f62d3096
# ╟─c898f2d0-b57b-4827-9a0f-43fd0dc0f7c9
# ╟─bbe9940c-1626-496b-80fa-00dfab30a031
# ╟─20cd6866-1508-4f00-956e-c64337207eca
# ╟─d7831b42-a57b-4f34-ae0a-61f00e9f9c1c
# ╟─6a2bb211-e55e-4973-bff2-310085340e70
# ╟─afd12f94-74fd-4f6f-a7b9-3aaf6a3e5ffe
# ╟─ca0483d3-fb9f-4035-a5dd-a6ab1da93605
# ╟─cc20226d-2d2e-4d5a-9936-d725e591eb45
