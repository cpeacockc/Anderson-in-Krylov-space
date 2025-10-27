using SparseArrays

function H_Anderson1D(L::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)
    
    #If an array of disorder is not given, one will be provided
    if isnothing(h)
        h = (W/2) .* (rand(L).*2 .-1) #random potentials
    end

    #Default is sparse format
    H_Anderson = spdiagm(-1 => (-t)*ones(L-1), 1 => (-t)*ones(L-1), 0=>h)

    #If periodic boundary conditions (default), connect the boundaries
    if periodic
    H_Anderson[1,L]=-t
    H_Anderson[L,1]=-t
    end
    return H_Anderson
end

function H_Anderson2D(L::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)

    H = spzeros(L^2,L^2)
    Id = spdiagm(ones(L))
        
        Big_Op = H_Anderson1D(L,0,t;periodic)
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,H_Anderson1D(L,0,t;periodic))
        H += Big_Op
        
    if isnothing(h)
        h = (W/2) .* (rand(L^2).*2 .-1) #random potentials
    end
    H += spdiagm(h)
    return H
end

function H_Anderson3D(L::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)

    H = spzeros(L^3,L^3)
    Id = spdiagm(ones(L))
        Big_Op = H_Anderson1D(L,0,t;periodic)
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,H_Anderson1D(L,0,t;periodic))
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,H_Anderson1D(L,0,t;periodic))
        H += Big_Op

    if isnothing(h)
        h = (W/2) .* (rand(L^3).*2 .-1) #random potentials
    end
    H += spdiagm(h)
    return H
end

function H_Anderson4D(L::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)

    H = spzeros(L^4,L^4)
    Id = spdiagm(ones(L))
    
        Big_Op = H_Anderson1D(L,0,t;periodic)
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,H_Anderson1D(L,0,t;periodic))
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,H_Anderson1D(L,0,t;periodic))
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,H_Anderson1D(L,0,t;periodic))
        H += Big_Op

    if isnothing(h)
        h = (W/2) .* (rand(L^3).*2 .-1) #random potentials
    end
    H += spdiagm(h)
    return H
end

function H_Anderson(d::Integer,L::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)

    if d==1
       return H_Anderson1D(L,W,t;h,periodic)
    elseif d==2
        return H_Anderson2D(L,W,t;h,periodic)
    elseif d==3
        return H_Anderson3D(L,W,t;h,periodic)
    elseif d==4
        return H_Anderson4D(L,W,t;h,periodic)
    else
        error("Incorrect dimension (only d=1-4 are currently included)")
    end
end