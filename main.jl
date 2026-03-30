using TensorKit, KrylovKit, Plots

include("tensors.jl")
include("plotting.jl")
include("tools.jl")
include("spectra.jl")

A = IsingTensor()

L = 12

untwisted(A, L)

psi_twisted(A, L)

L = 11
sigma_twisted(A, L)
