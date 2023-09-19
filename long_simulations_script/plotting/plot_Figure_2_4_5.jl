using Pkg
Pkg.activate(@__DIR__)
using DelimitedFiles, OrderedCollections
using StatsPlots, StatsBase
using LaTeXStrings
using RobustMeans
using Distributions

#! You might change that depending where your simulated data is stored
folder = "../cluster_script/g1sing_owen_data"

include("aux_plot_function.jl")
default(fontfamily="Computer Modern", linewidth=2, label=nothing)

# # Parameters
s = 0
m = 16 - s

base = 2

d = 10

N = base^m

n = 56 * (2^s)

M = 10 #10_000

ϵ = 0.05

δ = 3exp(-8) #  (4 or 8 in the exp)

# # Load points from cluster

QMC_method = [:MC, :NUS]
μ_QMC_unsort = point_from_folder(folder, N, d, ϵ, M, n, QMC_method)
μ_QMC = OrderedDict(name => μ_QMC_unsort[name] for name in unique([:MC;[keys(μ_QMC_unsort)...]])) # just making sure :MC is first

# One might need to resample Inf with the `resample_inf` as done in notebook `Figure_5.jl`.

# # Robust Estimators

pMN = 1 # Minsker Ndaoud hyperparameter

robust_names = ["Empirical Mean", "Huber", "Median of Mean", "Catoni", "Lee Valiant", "Minsker Ndaoud", "Trimmed Mean"]
estim() = [EmpiricalMean(), Huber(nothing), MedianOfMean(0), Catoni(nothing), LeeValiant(), MinskerNdaoud(pMN), TrimmedMean(0, 0)]
estimator_idx = 1:length(robust_names)-1
estimators = OrderedDict(short_name(robust_names[i]) => estim()[i] for i in estimator_idx)

# Long computation
results = OrderedDict{Symbol,OrderedDict}(name => computing_means(Float64.(reshape(U, n, M)), estimators, δ) for (name, U) in μ_QMC)

# Plot

# Figure 2 uses `true` to just display the quantile wisker. `false` displays the whole box + whisker
just_line = false

li = just_line ? "line" : "box"

seqs_name = length(μ_QMC) == 2 ? ["MC", "RQMC"] : collect(keys(μ_QMC))

sequence = keys(μ_QMC)

estimator = [:EM, :LV, setdiff(collect(keys(estimators)), [:EM, :LV])...] # Ensure it always starts with :EM and :LV

# ## Figure 2, 4 and 5

μ_exact = 1
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
            println("MSE(μ̂) = ", sqrt(var(val) + mean(val .- 1)^2))
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
    ylims!((0.999991e-4, 1.01e-1)) #! You might want to change that to adjust vertical axis
    yticks!(10.0 .^ (-9:-0))
    ylabel!(L"|\hat{\mu}_{N,n}-\mu|/\mu", tickfonthalign=:center)
end

# ## Figure 3

β = β_values[d]
const βₛ = β ./ [log(1 + s) for s in 1:d]
σ_exact = σ_F(βₛ)

qmc_seq = :NUS # :MC

begin
    fr = true # focus on central region ?
    framebox = fr ? "zoom" : "all"
    plot(thickness_scaling=2, size=(1000, 600))
    plot!(Normal(), label=L"\mathcal{N}(0,1)", c=:black, alpha=0.6)
    for (i, (name, s)) in enumerate(results[qmc_seq])
        W = √(n * N) * (s .- μ_exact) / σ_exact
        stephist!(W, alpha=0.6, norm=:pdf, label=string(name), c=i)
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