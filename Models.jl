using SparseArrays

function H_Anderson1D(L::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)
    
    #If an array of disorder is not given, one will be provided
    if isnothing(h)
        h = (W/2) .* (rand(L).*2 .-1) #random potentials
    end

    #Default is sparse format
    H_Anderson = spdiagm(-1 => (-t)*ones(L-1), 1 => (-t)*ones(L-1), 0=>h)

    #If periodic boundary conditions, connect the boundaries
    if periodic
    H_Anderson[1,L]=-t
    H_Anderson[L,1]=-t
    end
    return H_Anderson
end

function H_Anderson(d::Integer,L::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)
    
    #Total Hilbert space V
    V = L^d

    #fields given if not provided
    if isnothing(h)
        h = (W/2) .* (rand(V).*2 .-1) #random potentials
    end

    #Generate the 1d hopping model
    H_1D = H_Anderson1D(L,0,t,periodic) #periodic BCs by default

    H = spzeros(L^d,L^d)
    Id = spdiagm(ones(L))

    #Tensor product up to correct spatial dimension 
    for dim in 1:d
        Big_Ops = [Id for dim in 1:d]
        Big_Ops[d]=H_1D;Big_Op = Big_Ops[1]
        for j in 2:dim
            Big_Op = kron(Big_Op,Big_Ops[j])
        end
        H+=Big_Op
    end
    H += spdiagm(h)
    return H
end

