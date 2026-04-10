function IsingTensor()
    sp = Vect[IsingAnyon](:σ => 1)
    t = ones(ComplexF64, sp ⊗ sp ← sp ⊗ sp)
    block(t, IsingAnyon(:I)) .*= 1 + √(2)
    return t
end
function IsingBathroomTensor()
    sp = Vect[IsingAnyon](:σ => 1)
    t = ones(ComplexF64, sp ⊗ sp ← sp ⊗ sp)
    block(t, IsingAnyon(:I)) .*= 1 + √2
    block(t, IsingAnyon(:ψ)) .*= 1

    @planar T[-1 -2 -3 -4;-5 -6 -7 -8] := t[-1 1; -5 3] * t[-2 -3; 1 2] * t[3 4; -6 -7] * t[2 -4; 4 -8]

    F = isometry(fuse(sp, sp), sp ⊗ sp)

    @planar A[-1 -2;-3 -4] := F[-1;1 2] * F[-2;3 4] * T[1 2 3 4;5 6 7 8] * F'[5 6;-3] * F'[7 8;-4]
    return A / norm(A)
end
function TranslationMPO(sp::GradedSpace)
    id = complex(isomorphism(sp, sp))
    @planar t[-1 -2;-3 -4] := id[-1;-3] * id[-2;-4]
    return t
end
function half_braid(a::Sector, b::Sector, spv::GradedSpace)
    # Implements the halfbraid tensor for modular catys
    # spv = vertical space of transfer matrix

    theory = typeof(a) # Infer the anyon theory

    Z = Vect[theory](a => 1) ⊗ Vect[theory](b => 1)
    charge = fuse(Z)
    fusion = isometry(Z, charge)
    spvi = isomorphism(spv, spv)

    @plansor t[-1 -2;-3 -4] := spvi[-1;1] * τ[1 2;4 7] * conj(τ[5 6;7 3]) * spvi[6;-4] * fusion[4 5;-3] * conj(fusion[2 3;-2])
    return t
end
