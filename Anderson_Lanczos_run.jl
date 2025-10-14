using LinearAlgebra, Statistics,SparseArrays,PyPlot

#In diagonal basis we construct the Hamiltonian and Sz probe operator, then apply Lanczos

include("Lanczos_functions.jl");include("Models.jl")

Nsamples = 100 # number of samples to run
L = 5 # length of lattice
d=3 #spatial dimension
W=20 # disorder
t=1 # tunneling probability

Nsteps = 100 # Total number of steps of the Lanczos algorithm. In general, the Krylov space dimension is K_dim=V^2 -V + 1 where V=L^d
type_flag="dense-sparse" # options for lanczos operations are 'dense', 'sparse', or 'dense-sparse' which is recommended

# This code generates a probe which is the density on site L/2 (of each dimension, the exact origin)
site=Int(round(L/2)); Id = spdiagm(ones(L));
Diag = zeros(L)
Diag[site] = 1


bn_tot=zeros(Nsteps,Nsamples)

for k in 1:Nsamples
@show d,k, L, W, Nsteps
H=H_Anderson(d,L,W,t)


Probe = spdiagm(Diag)
for j in 2:d
    Probe = kron(Probe,Probe)
end

if type_flag=="dense"
    H = Matrix(H)
    Probe = Matrix(Probe)
elseif type_flag=="sparse"

elseif type_flag=="dense-sparse"
    Probe = Matrix(Probe)

else
    error("incorrect type_flag")
end

b = Lanczos(Probe,H,Nsteps)

b_tot[:,k]=b
end

b_avg=mean(b_tot,dims=2)[:]

ZeroMode_avg = ZeroMode_log_avg_from_b(b_tot)
ZeroMode_b_avg = ZeroMode_from_b(b_avg)
ZeroMode_diag = ZeroMode_from_diag(b_avg)

#plot the results
fig,axs=subplots()

ns=collect(1:length(ZeroMode_avg))
axs.plot(ns[1:2:end],ZeroMode_b_avg[1:2:end].^2,label="ZM from avg b")
axs.plot(ns[1:2:end],ZeroMode_avg[1:2:end].^2,label="ZM log-averaged over samples")
axs.plot(ns[1:2:end],ZeroMode_diag[1:2:end].^2,label="ZM from minimizing L^2")
axs.set_xlabel("Krylov space")
axs.set_ylabel("Zero mode probability density")
axs.set_yscale("log")
fig.savefig("Outs//Zero-modes.pdf")