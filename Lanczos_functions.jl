using LinearAlgebra, SparseArrays


#Liouvillian function for sparse matrices
function L_n(O::SparseMatrixCSC,H::SparseMatrixCSC,n::Integer)
    return H*O-O*H
end
#Liouvillian function for dense matrices
function L_n(O::AbstractMatrix,H::AbstractMatrix,n::Integer)
    prod = H*O
    return prod + ((-1)^(n+1))*adjoint(prod)
end
#Liouvillian function for dense O, sparse H (recommended)
function L_n(O::Matrix,H::SparseMatrixCSC,n::Integer) 
    prod = O*H
    prod += ((-1)^(n+1))*adjoint(prod)
    return -(prod)
end

#norm function for dense matrices
function Op_Norm(O::Matrix)
    return sqrt(tr(O'*O)/size(O)[1])
end
#norm function for sparse matrices
function Op_Norm(O::SparseMatrixCSC)
    return sqrt(sum(abs.(findnz(O)[3]).^2)/size(O)[1])
end
#Infinite temp inner product
function Op_Inner(A::AbstractArray,B::AbstractArray)
    return tr(A'*B)/size(A)[1]
end

function Lanczos(Probe::AbstractMatrix, H::AbstractMatrix, Nsteps::Integer)


    #Base Krylov vector 
    O = []

    #b1 is set to 1
    b = Float64[1]

    #If you want, track orthogonalization error
    #error_0 = Float64[]
    #error_1 = Float64[]
    
    #Define O0, normalize
    
    Probe/=Op_Norm(Probe)
    push!(O,Probe)
    LO_0 = L_n(O[1],H,2)

    #Define b2
    push!(b,Op_Norm(LO_0))
    
    #Define O1
    push!(O,LO_0/b[2])
    
    for n in 3:Nsteps
        #Construct new Krylov basis vector
        A_n = L_n(O[2],H,n) - b[n-1]*O[1]
        #Calculate Lanczos coefficient
        b_n = Op_Norm(A_n)
        push!(b,b_n)
        #Reset vectors, don't store entire basis in memory
        O[1]=O[2]
        O[2] = (A_n/b_n)

        #Optional during loop
        #push!(error_0,Op_Inner(Probe,O[2]))
        #push!(error_1,Op_Inner(LO_0/b[2],O[2]))
        #@show n,b,error_0, error_1
    end
    return b
end


#The following ZeroMode functions give the probability amplitudes of the zero modes in Krylov space
function ZeroMode_from_b(b)
    #This function assumes starting with b[1]=1 as in the previous code

    phi_n = zeros(length(b))
    phi_n[1]=1 #phi_0
    for n in 3:2:length(b)
        phi=1
        for i in 1:Int((n-1)/2)
            phi *=(b[2*i]/b[2*i+1]) 
        end
        phi_n[n]=phi
    end
    return phi_n
end

function ZeroMode_log_avg_from_b(b_tot)
    Nsamples = size(b_tot)[2]

    ZeroModes=zeros(length(ZeroMode_from_b(b_tot[:,1])),Nsamples)

    for j in 1:Nsamples
        ZeroModes[:,j]=ZeroMode_from_b(b_tot[:,j])
    end

    #We take the geometric mean
    ZeroMode_log_avgs=exp.(mean(log.(abs.(ZeroModes)),dims=2))
    return ZeroMode_log_avgs
end

#Construct Liouvillian in Krylov space
L_matrix(b) = Tridiagonal(b,zeros(length(b)+1),b)
