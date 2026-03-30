function plot_setup(L::Integer)
    title = "Transfer matrix L=$L"
    xaxis = "Topological spin"
    yaxis = "Scaling dimension"
    xticks = (-L / 2):1:(L / 2)
    yticks = [0, 1, 2, 3]
    legend = :bottomleft

    plot(
        title = title, xaxis = xaxis,
        xticks = xticks, yticks = yticks,
        yaxis = yaxis, legend = legend
    )

    ε = 0.1
    xlims!(-L / 2 - ε, L / 2 + ε)
    return ylims!(-ε, 2.9)
end
function plot_spectrum(spins::Vector{Float64}, dimensions::Vector{Float64}, sector::String; color = :red)
    ms = 4

    return scatter!(
        spins, dimensions,
        markershape = :diamond, markercolor = color,
        markersize = ms, label = sector
    )
end
function Ising_chars(sector::String)
    a, b = sector[1], sector[end]
    chars =
        Dict(
        '1' => [0.0, 2, 3, 3],
        'ψ' => [0.0, 1, 2, 3] .+ 1 / 2,
        'σ' => [0.0, 1, 2, 3, 3] .+ 1 / 16
    )
    χ, ϕ = chars[a], chars[b]
    s = [h₁ - h₂ for h₁ in χ for h₂ in ϕ]
    Δ = Float64[h₁ + h₂ for h₁ in χ for h₂ in ϕ]
    return s, Δ
end
function plot_analytical(sector; color = :red)
    s, Δ = Ising_chars(sector)

    ms, mal = 7, 0.35

    return scatter!(
        s, Δ,
        color = color, label = "$sector",
        markersize = ms, markeralpha = mal,
        markerstrokewidth = 0
    )
end
