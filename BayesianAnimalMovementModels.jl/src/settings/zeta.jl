
abstract type AbstractZeta end


struct ZetaDoNotUpdate <: AbstractZeta

end

#### DO UPDATE
struct ZetaDoUpdate <: AbstractZeta

    kmax::Int16

    emptyC::Vector{Int16}
    nonemptyC::Vector{Int16}

    n_obsC::Vector{Int16}

    n_emptyC::Vector{Int16}
    n_nonemptyC::Vector{Int16}

    zeta::Vector{Int16}

    n_itojC::Vector{Vector{Int16}}

    ZetaDoUpdate(kmax::Int16,emptyC::Vector{Int16},nonemptyC::Vector{Int16},n_obsC::Vector{Int16},n_emptyC::Vector{Int16},n_nonemptyC::Vector{Int16},zeta::Vector{Int16},n_itojC::Vector{Vector{Int16}}) = new(kmax,emptyC,nonemptyC,n_obsC,n_emptyC,n_nonemptyC,zeta,n_itojC)
end

function ZetaDoUpdate(zeta::Vector{Int16},kmax::Int16)

    if maximum(zeta)>kmax error("max of zeta must be less or equal than kmax") end
    if minimum(zeta)<1 error("min of zeta must be greater or equal to than 1") end

    emptyC      = Vector{Int16}(undef,kmax)
    nonemptyC   = Vector{Int16}(undef,kmax)

    n_obsC       = Vector{Int16}(undef,kmax)

    n_itojC     = Vector{Vector{Int16}}()

    for k in 1:kmax
        push!(n_itojC, zeros(Int16,kmax))
    end

    n_nonemptyC_1,n_emptyC_1 = Update_ObsInClust!(kmax,zeta,emptyC,nonemptyC,n_obsC,n_itojC)

    n_nonemptyC = Int16[n_nonemptyC_1]
    n_emptyC    =  Int16[n_emptyC_1]
    ZetaDoUpdate(kmax,emptyC,nonemptyC,n_obsC,n_emptyC,n_nonemptyC,zeta,n_itojC)

end

function Update_ObsInClust!(kmax::Integer,zeta::Vector{Int16},emptyC::Vector{Int16},nonemptyC::Vector{Int16},n_obsC::Vector{Int16},n_itojC::Vector{Vector{Int16}})::Tuple{Int16,Int16}

    for k in 1:kmax
        emptyC[k]       = Int16(0)
        nonemptyC[k]    = Int16(0)
        n_obsC[k]       = Int16(0)

        for k1 in 1:kmax
            n_itojC[k][k1] = Int16(0)
        end

    end

    zprev = 1
    for i in 1:(size(zeta,1)-1)
        n_itojC[zprev][zeta[i]] += one(zeta[i])
        n_obsC[zeta[i]]         += one(zeta[i])
        zprev                   = zeta[i]
    end

    iempty          = Int16(0)
    inonempty       = Int16(0)
    for k in 1:kmax
        if n_obsC[k] != 0
            inonempty += one(inonempty)
            nonemptyC[inonempty] = k
        else
            iempty += one(iempty)
            emptyC[iempty] = k
        end
    end
    return (inonempty , iempty)
end

function Update_ObsInClust!(zetastruct::ZetaDoUpdate)
    n_nonemptyC,n_emptyC = Update_ObsInClust!(zetastruct.kmax,zetastruct.zeta,zetastruct.emptyC,zetastruct.nonemptyC,zetastruct.n_obsC,zetastruct.n_itojC)

    zetastruct.n_emptyC[1,1]     = n_emptyC
    zetastruct.n_nonemptyC[1,1]  = n_nonemptyC
end
