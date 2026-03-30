function untwisted(A, L = 10)
    plot_setup(L)
    sector = [IsingAnyon(:I), IsingAnyon(:I)]
    spins, dimensions, _ = cft_spectrum(A, L, sector, amount = 12, Δ₀ = 0)
    plot_spectrum(spins, dimensions, "11")
    plot_analytical("11")

    sector = [IsingAnyon(:ψ), IsingAnyon(:ψ)]
    spins, dimensions, _ = cft_spectrum(A, L, sector, amount = 12, Δ₀ = 0)
    plot_spectrum(spins, dimensions, "ψψ", color = :green)
    return plot_analytical("ψψ", color = :green)
end
function psi_twisted(A, L = 10)
    L = 12
    plot_setup(L)

    sector = [IsingAnyon(:I), IsingAnyon(:ψ)]
    spins, dimensions, _ = cft_spectrum(A, L, sector, amount = 12, Δ₀ = 1 / 2)
    plot_spectrum(spins, dimensions, "1ψ")
    plot_analytical("1ψ")

    sector = [IsingAnyon(:ψ), IsingAnyon(:I)]
    spins, dimensions, _ = cft_spectrum(A, L, sector, amount = 12, Δ₀ = 1 / 2)
    plot_spectrum(spins, dimensions, "ψ1", color = :green)
    return plot_analytical("ψ1", color = :green)
end
function sigma_twisted(A, L)
    plot_setup(L)

    # sector = [IsingAnyon(:I), IsingAnyon(:σ)]
    # spins, dimensions, _ = cft_spectrum(A, L, sector, amount = 12, Δ₀ = 1/16)
    # plot_spectrum(spins, dimensions, "1σ")
    # plot_analytical("1σ")

    # sector = [IsingAnyon(:σ), IsingAnyon(:I)]
    # spins, dimensions, _ = cft_spectrum(A, L, sector, amount = 12, Δ₀ = 1/16)
    # plot_spectrum(spins, dimensions, "σ1", color=:green)
    # plot_analytical("σ1", color=:green)

    sector = [IsingAnyon(:ψ), IsingAnyon(:σ)]
    spins, dimensions, _ = cft_spectrum(A, L, sector, amount = 16, Δ₀ = 1 / 16)
    plot_spectrum(spins, dimensions, "ψσ")
    plot_analytical("ψσ")

    sector = [IsingAnyon(:σ), IsingAnyon(:ψ)]
    spins, dimensions, _ = cft_spectrum(A, L, sector, amount = 16, Δ₀ = 1 / 16)
    plot_spectrum(spins, dimensions, "σψ", color = :green)
    return plot_analytical("σψ", color = :green)
end
