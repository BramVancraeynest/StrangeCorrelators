using TensorKit, KrylovKit, Plots, ColorSchemes
using TNRKit

include("tensors.jl")
include("plotting.jl")
include("tools.jl")
include("spectra.jl")
include("anyonic_TNR.jl")

A = IsingTensor()

L = 12

untwisted(A, L)

psi_twisted(A, L)

L = 11
sigma_twisted(A, L)
