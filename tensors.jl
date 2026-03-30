function IsingTensor()
    sp = Vect[IsingAnyon](:σ => 1)
    t = ones(ComplexF64, sp ⊗ sp ← sp ⊗ sp)
    block(t, IsingAnyon(:I)) .*= 1 + √(2)
    return t
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
