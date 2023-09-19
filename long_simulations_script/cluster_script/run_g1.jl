# # Memory issues
# In Julia 1.9 on Linux (cluster) there is a known problem of garbage collection no working well. Hence long heavy computation can go "Out-Of-Memory". To avoid that one need to manually activate garbage collection. This was not a problem when I started experimenting on Cholesky cluster with Julia 1.7.
# https://discourse.julialang.org/t/poor-performance-of-garbage-collection-in-multi-threaded-application/75493/10?
# Not sure why 2^30 - 9
used_mem() = (println("$(round((Sys.total_memory()-Sys.free_memory())/2^30 -9))G used"))
# See also
# https://discourse.julialang.org/t/julia-killed-with-out-of-memory-error-on-linux-runs-fine-on-macos/102193/4

# # Settings

# ## Distributed
using Distributed
num_cores = 1 # Explicit declaration
## num_cores = parse(Int, ENV["SLURM_NPROCS"]) # on SLURM
## num_cores = parse(Int, ENV["NSLOTS"]) # on SGE

addprocs(num_cores)

@everywhere begin
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    Pkg.precompile()
end

println("Number of cores: ", nprocs()) #
println("Number of workers: ", nworkers())

# each worker gets its id, process id and hostname
for i in workers()
    id, pid, host = fetch(@spawnat i (myid(), getpid(), gethostname()))
    println(id, " ", pid, " ", host)
end
flush(stdout)

# ## Packages

@everywhere begin
    ## using SpecialFunctions, HypergeometricFunctions # used only for G_2sing function of the paper
    using QuasiMonteCarlo # add https://github.com/SciML/QuasiMonteCarlo.jl#master for latest versions
    using QMCGenerators: DigitalSeqB2G, Next
    using StatsBase: mean
end

using Random # for seeding even if it is useless for Distributed computing
Random.seed!(0)
using Dates

# ## Folders 

using DelimitedFiles

# # `Path` where you save your heavy computations
work_folder = "/home/dmetivier/work/g1sing_owen_data"
## work_folder = "g1sing_owen_data"

function create_folder(folder)
    if !isdir(folder)
        if !isdir(folder)
            mkdir(folder)
        else
            mkdir(folder)
        end
    end
end

# create_folder(folder)
create_folder(work_folder)

Today = today()

# ## Functions of the model
@everywhere begin
    F(x::AbstractArray, G, β) = prod(1 + β[i] * G(xᵢ) for (i, xᵢ) in enumerate(x))

    σ_F(β) = sqrt(prod(1 + βₖ^2 for βₖ in β) - 1)
end

# ## Parameters
@everywhere begin

    s = 0
    m = 16 - s

    base = 2

    d = 17

    N = base^m

    pad = 32

    n = 56 * (2^s)

    M = 10 #10_000

    const ϵ = 0.05

    β_values = [0.346574, 0.286243, 0.2615, 0.24686, 0.236718, 0.229051, 0.222928, 0.217849, 0.213519, 0.209752, 0.206421, 0.203437, 0.200736, 0.198269, 0.196, 0.1939, 0.191946, 0.190118, 0.188402, 0.186784]

    β = β_values[d] # could be codded as a function?
end

# Evaluate `const` of the conmputation
@everywhere begin
    const βₛ = β ./ [log(1 + s) for s in 1:d]

    G_1sing(x) = (sqrt(2ϵ) * (1 + 2ϵ)) / abs(1 - 2ϵ) * (x^(-1 / 2 + ϵ) - 2 / (1 + 2ϵ))

    f(x) = F(x, G_1sing, βₛ)

    operation(f, X) = mean(f(c) for c in eachcol(X))
end

# We choose our function such that μ = 1
μ_exact = 1

σ_exact = σ_F(βₛ)

# # Computation 

# Main loop
gc_rate = 0.1 # this value depends a lot on the size of the problem
# https://discourse.julialang.org/t/garbage-collection-not-aggressive-enough-on-slurm-cluster/61649/10?
function pmap_it(it, op)
    pmap(1:it.count) do i
        op(QuasiMonteCarlo.next!(it))
    end
end

function pmap_it(it, op, gc_rate)
    pmap(1:it.count) do i
        rand() < gc_rate && GC.gc()
        op(QuasiMonteCarlo.next!(it))
    end
end

# Dictionnary where we store MC and RQMC (NUS)
Upoints = Dict{Symbol,Vector}()
println("N=$(N) d=$(d) eps=$(round(ϵ, digits = 2)) M=$(M) n=$(n) pad=$(pad) $(Today)")
println("g2")

flush(stdout)

# ## MC

# ### Compute
used_mem()
@everywhere begin
    it_MC = DesignMatrix(N, d, RandomSample(), n * M, Float64)
end
@time Upoints[:MC] = pmap_it(it_MC, x -> operation(f, x))
used_mem()

# ### Save
for (key, val) in Upoints
    name = "U_N_$(N)_d_$(d)_eps_$(round(ϵ, digits = 2))_M_$(M)_n_$(n)_$(key).txt"
    writedlm(joinpath(work_folder, name), val)
end
flush(stdout)

# ## RQMC

# ### Design Matrix
@everywhere begin
    ## Add the method with Sobol starting at 0 from QMCGenerators.jl
    Base.@kwdef struct QMCGenDigitalSeqB2G <: QuasiMonteCarlo.DeterministicSamplingAlgorithm
        R::RandomizationMethod = NoRand()
    end

    function QuasiMonteCarlo.sample(n, d, ::QMCGenDigitalSeqB2G, T=Float64)
        return permutedims(Next(DigitalSeqB2G(d), n))
    end

    scrambling = OwenScramble(base=Int32(2), pad=Int32(32))
    key_scrambling = :NUS
    it_NUS = DesignMatrix(N, d, QMCGenDigitalSeqB2G(scrambling), n * M, Float64)
end

# ### The big computation
@time Upoints[key_scrambling] = pmap_it(it_NUS, x -> operation(f, x), gc_rate)
used_mem()
println(collect(keys(Upoints)))
flush(stdout)

# ### Save
for (key, val) in Upoints
    name = "U_N_$(N)_d_$(d)_eps_$(round(ϵ, digits = 2))_M_$(M)_n_$(n)_$(key).txt"
    writedlm(joinpath(work_folder, name), val)
end