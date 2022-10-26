module obj

using JuMP
using Ipopt

export calc_F,
       calc_F_λ_optimal_IPOpt,
       get_F_λ_opt

function calc_F(R,M)
    m = size(R)[1]
    f₁(w) = M'*w
    f₂(w) = sum([(R[t,:]'*w)^2 for t ∈ 1:m])/m
    f₃(w) = sum([(R[t,:]'*w)^3 for t ∈ 1:m])/m
    f₄(w) = sum([(R[t,:]'*w)^4 for t ∈ 1:m])/m
    return f₁, f₂, f₃, f₄
end

function calc_F_λ(R,M,λ)
    F = calc_F(R,M)
    F_λ(w) = -λ[1]*F[1](w) + λ[2]*F[2](w) - λ[3]*F[3](w) + λ[4]*F[4](w)
    return F_λ
end

function calc_F_λ_optimal_IPOpt(M,V,S,K,λ)
    n = length(M)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_cpu_time", 60.0)
    @variable(model, w[1:n] >= 0 )
    @constraint(model, w'*ones(n) == 1.0 )
    @NLobjective(model, Min, -λ[1]*sum(M[i]*w[i] for i in 1:n) +
                              λ[2]*sum(V[i,j]*w[i]*w[j] for i in 1:n for j in 1:n) -
                              λ[3]*sum(S[i,j,k]*w[i]*w[j]*w[k] for i in 1:n for j in 1:n for k in 1:n) + 
                              λ[4]*sum(K[i,j,k,l]*w[i]*w[j]*w[k]*w[l] for i in 1:n for j in 1:n for k in 1:n for l in 1:n) )
    optimize!(model)
    return model 
end

function get_F_λ_opt(M,V,S,K,λ)
    mod = calc_F_λ_optimal_IPOpt(M,V,S,K,λ)
    return F_λ_opt(extract_opt(mod)..., λ)
end

function extract_opt(mod)
    return value.(mod[:w]), 
           objective_value(mod),
           string(termination_status(mod)),
           solve_time(mod)
end

struct F_λ_opt
    optimizer 
    opt_val
    opt_stat
    sol_time
    λ
end

end
