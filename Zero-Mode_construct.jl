using LinearAlgebra, Statistics,SparseArrays
include("Lanczos_functions.jl");cd(@__DIR__)

#Define parameters
L = 5 # length of lattice
d=3 #spatial dimension
W=30 # disorder
t=1 # tunneling probability
Nsteps=200

#Load in some data
#path_out="C://MyDrive//Documents//A-Physics-PhD//Dries-Research//Code//Anderson_krylov//Outs//Lanczos//"*H_flag
#cd(path_out)
cd("Outs//Lanczos_data//b-d$d-L$L-W$W-t$t")
Nsamples = length(readdir()) #Take all samples in the folder
b_tot = zeros(Nsteps,Nsamples)
for k in 1:Nsamples
    file="b$k.txt"
    for line in eachline(file)
    eval(Meta.parse(line))
    end
    b_tot[:,k]=b 
end

#Calculate the average
b_avg=mean(b_tot,dims=2)[:]

#Construct the zero modes
#First from each sample, averaging after
ZeroMode_avg = ZeroMode_log_avg_from_b(b_tot)
#Next from average b
ZeroMode_b_avg = ZeroMode_from_b(b_avg)

#plot the results, comparing the two zero modes and checking
#
using PyPlot
fig,axs=subplots()

ns=collect(1:length(ZeroMode_avg)) .^(1/(2*d)) #Rescale as predicted
axs.plot(ns[1:2:end] ,ZeroMode_b_avg[1:2:end].^2,label="ZM from average b")
axs.plot(ns[1:2:end],ZeroMode_avg[1:2:end].^2,label="ZM log-averaged over samples")
axs.set_xlabel("Krylov space n^(1/$(2d))")
axs.set_ylabel("Zero mode probability density")
axs.set_yscale("log")
fig
fig.savefig("..//Outs//Zero-modes.pdf")
