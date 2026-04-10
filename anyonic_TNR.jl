# finalize with half-braid
function custom_finalize!(scheme::TNRKit.simple_scheme)
    return scheme.T / norm(scheme.T)
end

num_steps = 8
num_dims = 12
start_iter = 3 # Compute data starting at this RG step
L = 2 # Circumference of transfer matrix final ED
rnk = 16 # truncrank

t = IsingBathroomTensor()

scheme = TRG(t)
finalizer = Finalizer(custom_finalize!, AbstractTensorMap)
data = run!(scheme, truncrank(rnk), maxiter(20), finalizer; finalize_beginning = true, verbosity = 1); # works on branch

# Post processing
sector = [IsingAnyon(:σ), IsingAnyon(:σ)]

dimensions_flow = zeros(num_dims, num_steps)
spins = similar(dimensions_flow)

for index in 1:num_steps
    step = start_iter + index
    T_step = data[step]
    @info "Computating data for RG step $(step)"
    ss, dimensions, _ = cft_spectrum(T_step, L, sector, amount = num_dims, Δ₀ = 0, mom = false)
    dimensions_flow[:, index] = dimensions[1:num_dims]
    spins[:, index] = ss[1:num_dims]
end


plot_setup_RGflow(L, start_iter, rnk)

for i in 1:num_dims
    plot!(1:num_steps, dimensions_flow[i, :], marker = :x, label = "$i", palette = :acton25)
end
plot!()
