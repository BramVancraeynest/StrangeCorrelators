function action(A::AbstractTensorMap{S₁, N₁, N₂}, halfbraid::AbstractTensorMap{S₃, N₃, N₄}, x::AbstractTensorMap{S₂, d₁, d₂}) where {S₁, N₁, N₂, S₃, N₃, N₄, S₂, d₁, d₂}
    transl = TranslationMPO(codomain(A)[1])

    t = apply_mpo(A, halfbraid, x)
    t = apply_mpo(A, halfbraid, t)

    t = apply_mpo(transl, halfbraid, t)

    return t
end

tensorexpr(name::Symbol, indout, indin) = Expr(:typed_vcat, name, Expr(:row, indout...), Expr(:row, indin...))
@generated function apply_mpo(A::AbstractTensorMap{T₁, S₁, N₁, N₂}, halfbraid::AbstractTensorMap{T₂, S₂, N₃, N₄}, x::AbstractTensorMap{T₃, S₃, N, 1}) where {T₁, S₁, N₁, T₂, N₂, S₂, N₃, N₄, T₃, S₃, N}
    N > 0 || error("undefined behaviour for length 0")

    out_part = tensorexpr(:y, -1:-1:-N, (-N - 1))

    in_part = Expr(:call, :*, tensorexpr(:x, 1:2:(2N - 1), 2N + 2), tensorexpr(:A, (-1, 2), (2N + 1, 1)))
    for i in 2:N
        push!(in_part.args, tensorexpr(:A, (-i, 2i), (2i - 2, 2i - 1)))
    end

    push!(in_part.args, tensorexpr(:halfbraid, (2N + 1, 2N + 2), (-N - 1, 2N)))

    return :(@plansor $out_part := $in_part)
end
function cft_spectrum(A::AbstractTensorMap{S, N₁, N₂}, L::Integer, sector::Vector; Δ₀ = 0::Number, amount = 10::Integer, n_vec = 1::Integer) where {S, N₁, N₂}
    sp = domain(A)[2] # Infer spaces from transfer matrix tensor
    spv = codomain(A)[2]

    charge = Vect[typeof(sector[1])](x => 1 for x in sector[1] ⊗ sector[2])
    init = rand(ComplexF64, sp^L, charge)

    norm(init) > 1.0e-12 || error("Initial vector has norm 0, length $L not compatible with charge $charge !")

    hb = half_braid(sector..., spv)

    spectrum, vecs, _ = eigsolve(x -> action(A, hb, x), init, krylovdim = 20, amount, :LM)

    spectrum ./= norm(spectrum[1]) # Get rid of non-universal contributions
    spectrum .*= exp(-2π / L * 2Δ₀) # Need double Δ₀
    spins = angle.(spectrum) .* L / (2 * π)
    dimensions = -L / (4 * π) .* log.(norm.(spectrum))

    return spins, dimensions, vecs
end
