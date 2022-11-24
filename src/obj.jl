module obj

using JuMP
using Ipopt

using Random
using Combinatorics
using Distributed

export get_F_λ_opt,
       calc_F_λ_optimal_IPOpt,
       calc_F_λ_optimal_sparse_IPOpt
 

struct F_λ_opt
    domain
    optimizer 
    opt_val
    opt_stat
    sol_time
    λ
end

function get_F_λ_opt(M,V,S,K,λ; box_size=0, sub=0)
        domain = box_size == 0 ? "simplex" : "box"
        if sub == 0
            mod = calc_F_λ_optimal_IPOpt(M,V,S,K,λ,box_size)
            return F_λ_opt(domain, extract_opt(mod)..., λ)
        else
            mod, sup, comp_time = calc_F_λ_optimal_sparse_IPOpt(M,V,S,K,λ,sub, box_size)
            n   = length(M)
            return F_λ_opt(domain, extract_opt(mod,sup,n)..., comp_time,λ)
        end
end

#---------------------------------------------------------------------------------------------------------
function calc_F_λ_optimal_sparse_IPOpt(M,V,S,K,λ,k; box_size=0)
    n = length(M)
    subs = collect(powerset([1:n...],k,k))

    mod_vec = []
    @distributed for s in subs
        push!(mod_vec, calc_F_λ_optimal_IPOpt(M[s],V[s,s],S[s,s,s],K[s,s,s,s],λ,box_size))
    end

    opt_vec   = map(mod -> JuMP.objective_value(mod), mod_vec)
    sp        = sortperm(opt_vec)[end]
    comp_time =  sum(map(mod -> JuMP.solve_time(mod), mod_vec))
    return mod_vec[sp[1]], subs[sp[1]], comp_time
end

function calc_F_λ_optimal_IPOpt(M,V,S,K,λ,box_size=0)
    n = length(M)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_cpu_time", 60.0) # cut-off time

    if box_size > 0 # box constraints
        @variable(model, w[1:n]   >= -box_size )
        @constraint(model, w[1:n] .<=  box_size)
    else # simplex constraints
        @variable(model, w[1:n] >= 0 )                      
        @constraint(model, w'*ones(n) == 1.0 )      
    end        

    @NLobjective(model, Min, -λ[1]*sum(M[i]*w[i] for i in 1:n) +
                              λ[2]*sum(V[i,j]*w[i]*w[j] for i in 1:n for j in 1:n) -
                              λ[3]*sum(S[i,j,k]*w[i]*w[j]*w[k] for i in 1:n for j in 1:n for k in 1:n) + 
                              λ[4]*sum(K[i,j,k,l]*w[i]*w[j]*w[k]*w[l] for i in 1:n for j in 1:n for k in 1:n for l in 1:n) )
    optimize!(model)
    return model 
end

#---------------------------------------------------------------------------------------------------------
function extract_opt(mod)
    return value.(mod[:w]), 
           objective_value(mod),
           string(termination_status(mod)),
           solve_time(mod)
end

function extract_opt(mod,sub,n)
    w = zeros(n)
    w[sub] = value.(mod[:w])
    return w, 
           objective_value(mod),
           string(termination_status(mod))
end

function calc_F(R,M)
    m = size(R)[1]
    f₁(w) = M'*w
    f₂(w) = sum([(R[t,:]'*w)^2 for t ∈ 1:m])/m
    f₃(w) = sum([(R[t,:]'*w)^3 for t ∈ 1:m])/m
    f₄(w) = sum([(R[t,:]'*w)^4 for t ∈ 1:m])/m
    return f₁, f₂, f₃, f₄
end

end
