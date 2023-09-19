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
m_max = 16 # 13

d = 10

M = 10 #10_000

ϵ = 0.05

δ = 3exp(-8) #  (4 or 8 in the exp)

base = 2

# # Load points from cluster

QMC_method = [:NUS]

# One might need to resample Inf with the `resample_inf` as done in notebook `Figure_5.jl`.

# # Robust Estimators

# In the paper we show that plot only for :NUS with Lee Valiant.
estimator = LeeValiant()
name_esti = :LV

sequence = [:NUS]

seqs_name = [string.("RQMC", " + Robust")]

estimators = Dict(name_esti => estimator)

# `Number_of_files` is the number of constant budget B = n*N you want to test. In the paper we do `Number_of_files = 4` 
Number_of_files = 2

μ_exact = 1

# # Figure 6

begin
    x_coord = zeros(Int, 0)
    y_coord = zeros(0)
    z_coord = zeros(Int, 0)

    plot(thickness_scaling=2)
    expo = 0:Number_of_files-1
    for (id, ex) in enumerate(expo)
        m = m_max - ex
        N = base^m
        n = 56 * (2^ex)
        Upoints = point_from_folder(folder, N, d, ϵ, M, n, sequence)

        ## resample_inf!(Upoints, DesignMatrix(N, d, SobolSample(OwenScramble(base=base, pad=pad)), 1), x -> operation(F_test, x))
        
        results = OrderedDict{Symbol,OrderedDict}(name => computing_means(reshape(U, n, M), estimators, δ) for (name, U) in Upoints)

        for (seq_num, s) in enumerate(sequence)
            for (kk, estimate) in enumerate(keys(estimators))
                kk = ex + 1

                val = results[s][estimate] / μ_exact
                ΔM = abs.(val .- 1)
                append!(y_coord, ΔM)
                append!(x_coord, fill(kk, length(ΔM)))
                append!(z_coord, fill(seq_num, length(ΔM)))
            end
        end
    end
    groupedboxplot!(x_coord, y_coord, group=z_coord, label=permutedims(seqs_name), size=(1200, 600), alpha=1, c=2, ms=1.8,
        whisker_range="$(round(1-δ, digits = 6))"
    )
    estimator_name = map(expo) do ex
        m = m_max - ex
        N = base^m
        n = 56 * (2^ex)
        return string(latexstring("N = "), N, "\n", latexstring("n = "), n)
    end
    xticks!(1:length(estimator_name), string.(estimator_name), background_color_legend=nothing)
    yaxis!(:log10, yminorticks=9, minorgrid=:y, legend=(0.35, 0.93), minorgridlinewidth=1.2)
    ylims!((0.999991e-5, 2.01e-3))
    yticks!(10.0 .^ (-9:-0))
    ylabel!(L"|\hat{\mu}_{N,n}-\mu|/\mu", tickfonthalign=:center)
end

# ## Save

# Save the figure and crop it to avoid white spaces appearing in gr()
li = today()

begin # save the figure and crop it to avoid white spaces appearing in gr()
    save_name = "save/boxplot_N_vs_n_d_$(d)_esp_$(round(ϵ, digits = 2))_delta_$(round(δ, digits = 3))_$(string(estimator...))_M_$(M)_$(li)"
    save_name = replace(save_name, "." => "p")
    savefig(string(save_name, ".pdf"))
    run(`pdfcrop $(string(save_name,".pdf"))`) # Petit délire pour croper proprement la figure 
    mv(string(save_name, "-crop", ".pdf"), string(save_name, ".pdf"), force=true)
end