module obj

using JuMP
using Ipopt


# using Random
# using Combinatorics

include("spaces.jl")                ; using .spaces

export get_F_λ_opt,
       calc_F_λ_optimal_IPOpt,
       calc_F_λ_optimal_FISTA,
    #    calc_F_λ_optimal_sparse_IPOpt,
       calc_F_λ,
       calc_∇F_λ
 
struct F_λ_opt
    domain
    optimizer 
    opt_val
    opt_stat
    sol_time
    λ
end

# function get_F_λ_opt(M,V,S,K,λ, box_size, sub)
#         domain = box_size == 0 ? "simplex" : "box"
#         # if sub == 0
#             mod = calc_F_λ_optimal_IPOpt(M,V,S,K,λ,box_size)
#             return F_λ_opt(domain, extract_opt(mod)..., λ)
#         # else
#             # mod, sup, comp_time = calc_F_λ_optimal_sparse_IPOpt(M,V,S,K,λ,sub, ; box_size=box_size)
#             # n   = length(M)
#             # return F_λ_opt(domain, extract_opt(mod,sup,n)..., comp_time,λ)
#         # end
# end


function get_F_λ_opt(pd,λ;x₀,iter=0, box_size=0, sub=0)
    domain = dom_choose(box_size,sub)
    if iter > 0
        t = @timed (x, fx) = calc_F_λ_optimal_FISTA(pd.R,pd.M,pd.V,pd.S_mat,pd.K_mat,λ,x₀,iter,box_size)
        F_λ_opt(domain,x,fx,"Unknown",t.time,λ)
    else
        mod = calc_F_λ_optimal_IPOpt(pd.M,pd.V,pd.S,pd.K,λ,box_size)
        F_λ_opt(domain, extract_opt(mod)..., λ)
    end
end
dom_choose(box_size=0, sub=0) = (box_size == 0 ? "simplex" : "box")*"_supp_"*(sub == 0 ? "dense" : string(sub))

function calc_F_λ_optimal_IPOpt(M,V,S,K,λ,box_size)
    n = length(M)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_cpu_time", 60.0) # cut-off time
    set_optimizer_attribute(model, "warm_start_init_point_", true)
    # set_optimizer_attribute(model, "warm_start_initializer", true)
    

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

function calc_F_λ_optimal_FISTA(R,M,V,S_mat,K_mat,λ,x₀,iter,box_size)   
    f  = obj.calc_F_λ(M, V, S_mat, K_mat, λ) 
    ∇f = obj.calc_∇F_λ(R, M, V, λ)
    x = FISTA(∇f,x₀,iter,box_size)
    fx = f(x)
    return x, fx 
end
function FISTA(∇f,x₀,iter,box_size)
    # initiallize
    yₖ = xₖ = x₀
    tₖ = 1
    L = Lₖ = 4

    if box_size == 0
        proj = spaces.Euc_proj_simp
    else
        proj = x -> spaces.Euc_proj_box(x,box_size)
    end

    # iterate 
    for _ = 1:iter
        # Lₖ  = update_Lₖ(f, ∇f, yₖ, Lₖ)
        tₖ₊₁ = (1+sqrt(1+4*tₖ^2))/2 

        # xₖ₊₁ = prox((1/Lₖ), τ, yₖ - (1/Lₖ)*∇f(yₖ))
        xₖ₊₁ = proj(yₖ - (1/Lₖ)*∇f(yₖ))
        yₖ₊₁ = xₖ₊₁ + ((tₖ-1)/(tₖ₊₁))*(xₖ₊₁-xₖ)
        
        yₖ =yₖ₊₁
        xₖ = xₖ₊₁
        tₖ = tₖ
    end
    return xₖ
end
function update_Lₖ(f,∇f, yₖ, Lₖ)
    η = 2   # η > 1
    iₖ = 1

    Tyₖ = (1/(Lₖ*η^(iₖ)))*spaces.Euc_proj_simp(yₖ) 
    Tyₖ_min = Tyₖ - yₖ
    
    while  f(Tyₖ) ≤ f(yₖ) + ∇f(yₖ)'*(Tyₖ_min) + (Lₖ/2)*sqrt(Tyₖ_min'*Tyₖ_min) 
        Lₖ = Lₖ*η
        Tyₖ = (1/(Lₖ*η^(iₖ)))*spaces.Euc_proj_simp(yₖ) 
        Tyₖ_min = Tyₖ - yₖ 
        iₖ+= 1
    end

    return Lₖ
end

#---------------------------------------------------------------------------------------------------------
# function calc_F_λ_optimal_sparse_IPOpt(M,V,S,K,λ,k; box_size=0)
#     n = length(M)
#     subs = collect(powerset([1:n...],k,k))

#     mod_vec = []
#     @distributed for s in subs
#         push!(mod_vec, calc_F_λ_optimal_IPOpt(M[s],V[s,s],S[s,s,s],K[s,s,s,s],λ,box_size))
#     end

#     opt_vec   = map(mod -> JuMP.objective_value(mod), mod_vec)
#     sp        = sortperm(opt_vec)[end]
#     comp_time =  sum(map(mod -> JuMP.solve_time(mod), mod_vec))
#     return mod_vec[sp[1]], subs[sp[1]], comp_time
# end
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

function calc_F_λ(M, V, S_mat, K_mat,λ)
    f₁, f₂, f₃, f₄ = calc_F(M, V, S_mat, K_mat)
    return w ->  [-λ[1],λ[2],-λ[3],λ[4]]'*[f₁(w), f₂(w), f₃(w), f₄(w)]
end

function calc_F(M, V, S_mat, K_mat)
    f₁(w) = M'*w
    f₂(w) = w'*V*w
    f₃(w) = kron(w,w)'*S_mat*w
    f₄(w) = kron(w,w)'*K_mat*kron(w,w)
    return f₁, f₂, f₃, f₄
end

function calc_∇F_λ(R,M,V,λ)
    ∇f₁, ∇f₂, ∇f₃, ∇f₄ = calc_∇F(R,M,V)
    return w ->  [-λ[1],λ[2],-λ[3],λ[4]]'*[∇f₁(w), ∇f₂(w), ∇f₃(w), ∇f₄(w)]
end

function calc_∇F(R,M,V)
    m = size(R)[1]
    ∇f₁(w) = M
    ∇f₂(w) = 2*V*w 
    ∇f₃(w) = 3*sum([(R[t,:]-M)*((R[t,:]-M)'*w)^2 for t ∈ 1:m])/m
    ∇f₄(w) = 4*sum([(R[t,:]-M)*((R[t,:]-M)'*w)^3 for t ∈ 1:m])/m
    return ∇f₁, ∇f₂, ∇f₃, ∇f₄
end


end






# function calc_F(R,M)
#     m = size(R)[1]
#     f₁(w) = M'*w
#     f₂(w) = sum([(R[t,:]'*w)^2 for t ∈ 1:m])/(m-1)
#     f₃(w) = sum([(R[t,:]'*w)^3 for t ∈ 1:m])/m
#     f₄(w) = sum([(R[t,:]'*w)^4 for t ∈ 1:m])/m
#     return f₁, f₂, f₃, f₄
# end
# function calc_F_λ(R,M,λ)
#     f₁, f₂, f₃, f₄ = calc_F(R,M)
#     return w ->  [-λ[1],λ[2],-λ[3],λ[4]]'*[f₁(w), f₂(w), f₃(w), f₄(w)]
# end