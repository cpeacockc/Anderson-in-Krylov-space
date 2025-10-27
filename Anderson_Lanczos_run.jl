using LinearAlgebra, Statistics,SparseArrays
cd(@__DIR__)
#In diagonal basis we construct the Hamiltonian and probe operator (Z on a single site), then apply Lanczos

include("Lanczos_functions.jl");include("Models.jl")

L = 10 # length of lattice
d=3 #spatial dimension
W=40 # disorder
t=1 # tunneling probability
Nsteps = 1000 # Total number of steps of the Lanczos algorithm. In general, the Krylov space dimension is K_dim=V^2 -V + 1 where V=L^d
Nsamples = 100
type_flag="dense-sparse" # options for lanczos operations are 'dense', 'sparse', or 'dense-sparse' which is recommended

#Create output path for storing coefficient samples
mkpath("Outs//Lanczos_data//b-d$d-L$L-W$W-t$t")
cd("Outs//Lanczos_data//b-d$d-L$L-W$W-t$t")


#k=1 #Tracks the realization/sample number if running single samples
for k in 1:Nsamples

    #Generates Anderson Hamiltonian with dimension d and random fields, default periodic BCs
    #H=H_Anderson(d,L,W,t;h=h;periodic=false)
    
    H=H_Anderson(d,L,W,t)

    #Generates a probe which is the density on site L/2 (of each dimension, the exact origin)
    site=Int(round(L/2)); Id = spdiagm(ones(L));
    Diag = zeros(L)
    Diag[site] = 1 
    Probe = spdiagm(Diag)
    
    for j in 2:d
        Probe = kron(Probe,spdiagm(Diag))
    end

    #Sets the storage type of the Probe and H
    if type_flag=="dense" #In this both the probe and H are stored as dense.
        H = Matrix(H)
        Probe = Matrix(Probe)

    elseif type_flag=="sparse" #In this case both probe and H are stored as sparse. Scales worse for memory, better for speed

    elseif type_flag=="dense-sparse" #Default: Probe is stored as dense but H is stored as sparse 
        
        Probe = Matrix(Probe)
    else
        error("incorrect type_flag")
    end

    #Generate the Lanczos coefficients
    b = Lanczos(Probe,H,Nsteps)

    #Store the coefficients in a text file
    file = "b$k.txt"
    f = open(file,"w")
    println(f,"b = $b")
    close(f)
end
