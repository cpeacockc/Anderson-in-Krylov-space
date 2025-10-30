# Anderson Localization in Krylov space
Julia code for calculating Lanczos coefficients for the Anderson model in Krylov space, as well as some functionality to construct the integrals of motion, which are the zero modes of the Liouvillian. *Hamiltonians and operators are written in single particle basis*

### Models.jl 
Stores the Hamiltonians for the various anderson models

### Lanczos_functions.jl 
Stores all other functions needed

### Anderson_Lanczos_run.jl 
Example code calculating coefficients in the Anderson model, running for a certain number of samples

### Zero-Mode_construct.jl 
Example code of constructing zero-modes and plotting. Example plotting code is written for PyPlot.jl but is default commented out

## Quickstart Functions

### Create Anderson Hamiltonian of dimension d, lattice length N, disorder W, tunneling probability t
```julia
d=3;N=5;t=1;W=30
H = H_Anderson(d,N,t,w)
```

### Create the probe operator (here a Pauli Z operator in the middle of the torus)
```julia
site=Int(round(L/2)); Id = spdiagm(ones(L));
Diag = zeros(L)
Diag[site] = 1 
Probe = spdiagm(Diag)
for j in 2:d
  Probe = kron(Probe,spdiagm(Diag))
end
```

### Run the Lanczos algorithm to calculate the Lanczos coefficients for this H and Probe
Storing the Hamiltonian as sparse and the probe as dense is generally recommended
```julia
Nsteps=100
b = Lanczos(Matrix(Probe),H,Nsteps)
```

### Calculate the average zero-modes, or conserved quantities
After collecting some samples of coefficients, stored in ``b_tot[Nsteps,Nsamples]'', we can check the form of the average zero-modes. 
First average the coefficients first over samples and then calculate the zero mode
```julia
b_avg=mean(b_tot,dims=2)[:]
ZeroMode_b_avg = ZeroMode_from_b(b_avg)
```
We can then compare this to the zero modes calculated for each sample, then log-averaged
```julia
ZeroMode_avg = ZeroMode_log_avg_from_b(b_tot)
```

### Plot the zero-modes
Example of plotting both zero-modes to compare, and plotting against the predicted scaling in the localized regime: ZeroMode(n) ~ exp(-n^(1\2d))
```julia
using PyPlot
fig,axs=subplots()

ns=collect(1:length(ZeroMode_avg)) .^(1/(2*d)) #Rescale as predicted
axs.plot(ns[1:2:end] ,ZeroMode_b_avg[1:2:end].^2,label="ZM from average b")
axs.plot(ns[1:2:end],ZeroMode_avg[1:2:end].^2,label="ZM log-averaged over samples")
axs.set_xlabel("Krylov space n^(1/$(2d))")
axs.set_ylabel("Zero mode probability density")
axs.set_yscale("log")
fig
```
Note that the zero mode is exactly zero on every other site
