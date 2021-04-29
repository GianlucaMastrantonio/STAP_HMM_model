
function compute_AppCondMeanAndVariance_MvNormal(IndexA::Vector{T},CovInv::Matrix{Float64} ) where{T<:Integer}

    nvar        = size(CovInv,1)
    IndexB      = T[1:nvar;]
    IndexB      = deleteat!(IndexB,IndexA)

    CovInv_AgivenB = CovInv[IndexA,IndexA]
    CovAgivenB     = inv(CovInv_AgivenB)

    CovABInvCovB      = -CovAgivenB*view(CovInv,IndexA,IndexB)

    return (CovABInvCovB, CovAgivenB ,CovInv_AgivenB )


end


function compute_CondMeanAndVariance_MvNormal(IndexA::Vector{T},Mean::Vector{Float64}, CovInv::Matrix{Float64}, Obs::Vector{Float64} ) where {T<:Integer}

    nvar        = size(CovInv,1)
    IndexB      = T[1:nvar;]
    IndexB      = deleteat!(IndexB,IndexA)

    CovABInvCovB, CovAgivenB ,CovInv_AgivenB = compute_AppCondMeanAndVariance_MvNormal(IndexA,CovInv )
    MuAgivenB      = Mean[IndexA] + CovABInvCovB*(view(Obs,IndexB)- view(Mean,IndexB) )

    return ( MuAgivenB, CovAgivenB ,CovInv_AgivenB )


end


function logsumexp(X::Vector{Float64})
    alpha = -Inf
    r = 0.0
    for x = X
        if x <= alpha
            r += exp(x - alpha)
        else
            r *= exp(alpha - x)
            r += 1.0
            alpha = x
        end
    end
    return log(r) + alpha
end
function sample_discretevar(ProbNonNorm::Vector{Float64})
#X = ProbNonNorm
    App = exp.(ProbNonNorm .-logsumexp(ProbNonNorm))
    #return findall(rand(Multinomial(1,App/sum(App) )) .==1)[1];
    return findall(cumsum(App).>=rand(Uniform(0.0,sum(App))))[1]
    #return 1

end
#
# int Class_Utils::sample_DiscreteVar(double *logprob_NonNormalized, int K)
# {
#     int k;
#     double sum = 0.0;;
#     double prob[K];
#     Class_Utils::log_sum_exp(&sum, K, logprob_NonNormalized);
#
#     for(k=0;k<K;k++)
#     {
#         prob[k] = exp(logprob_NonNormalized[k]-sum);
#     }
#
#     double u = runif(0.0,1.0);
#     //REprintf("u %f ",u);
#     sum = 0.0;
#     k = -1;
#     do{
#         k++;
#         sum += prob[k];
#         //REprintf("%f ",sum );
#     }while(sum<u);
#
#     //REprintf("K=%i\n",k);
#     return(k);
# }
