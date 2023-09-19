function short_name(name::AbstractString)
    s = split(name)
    return length(s) > 1 ? Symbol(first.(s)...) : Symbol(uppercase(string(first.(s, 2)...)))
end

# # Load points from cluster

searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))

"""
    point_from_folder(folder, N, d, ϵ, M, n, QMC_method)
Find all files in folder matching the parameters
"""
function point_from_folder(folder, N, d, ϵ, M, n, QMC_method)
    file_name = "U_N_$(N)_d_$(d)_eps_$(round(ϵ, digits = 2))_M_$(M)_n_$(n)"
    file_names = searchdir(folder, file_name)
    idx_file = [findfirst(true .== [occursin(string(QMC), file) for file in file_names]) for QMC in QMC_method]
    return Dict(key => vec(readdlm(joinpath(folder, file_names[idx_file[i]]))) for (i, key) in enumerate(QMC_method))
end

# # Robust

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

function RobustMeans.mean(A::AbstractArray, δ::Real, Estimator::Catoni{<:Nothing}, kwargs...)
    n = length(A)
    α = RobustMeans.α_Catoni(δ, n, std(A))
    z = RobustMeans.Z_Estimator(α, RobustMeans.ψ_Catoni)
    return mean(A, z; kwargs...)
end

function RobustMeans.mean(A::AbstractArray, δ::Real, Estimator::Huber{<:Nothing}, kwargs...)
    n = length(A)
    α = RobustMeans.α_Huber(δ, n, std(A))
    z = RobustMeans.Z_Estimator(α, RobustMeans.ψ_Huber)
    return mean(A, z; kwargs...)
end

# # Plot 

import StatsPlots: notch_width, _cycle

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

# # `F` function std σ_exact

σ_F(β) = sqrt(prod(1 + βₖ^2 for βₖ in β) - 1)

β_values = [0.346574, 0.286243, 0.2615, 0.24686, 0.236718, 0.229051, 0.222928, 0.217849, 0.213519, 0.209752, 0.206421, 0.203437, 0.200736, 0.198269, 0.196, 0.1939, 0.191946, 0.190118, 0.188402, 0.186784]

