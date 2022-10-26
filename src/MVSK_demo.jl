include("stock_data.jl")            ; using .stock_data          
include("λ_char.jl")                ; using .λ_char 
include("pgd.jl")                   ; using .pgd       
include("spaces.jl")                ; using .spaces
include("obj.jl")                   ; using .obj

include(pwd()*"/test/runtests.jl") 

#-------------------------------------------- Data -------------------------------------------- 
pd = stock_data.load_processed_data()
pd_5 = stock_data.load_processed_data([1:5...],[1:5...])
# typeof(pd)
# fieldnames(typeof(pd))
#-------------------------------------------- λ₃ characterizations -------------------------------------------- 
λ₂, λ₃, λ₄ = spaces.gen_simp_ele()[1:3]

λ₃_bound_all  = λ_char.λ₃_bound(λ₂, λ₄)
λ₃_bound_simp = λ_char.λ₃_bound_Δ(λ₂, λ₄, pd)

# These all take long and give bad estimates.
λ₃_bound_I    = λ_char.λ₃_bound_IPOpt(λ₂, λ₄, pd_5)
λ₃_bound_s    = λ_char.λ₃_bound_SOS(λ₂, λ₄, pd_5, 4) #λ₃ SOS bounds (NEEDS TO BE TESTED)# One test is if it corresponds to λ₃_bound_all  when not restricted to the simplex
λ₃_bound_m    = λ_char.λ₃_bound_MOM(λ₂, λ₄, pd_5, 2)  # λ₃ moment bounds, very bad bounds
# λ₃_bound_g  = λ_char.λ₃_bound_pgd(λ₂, λ₄, pd.R[:,1:6],2)   # has not been coded before
#-------------------------------------------- Objective values-------------------------------------------- 
pd  = stock_data.get_proc_data()
#----------------------------------------------------------------------------------------Single point
λ = spaces.gen_simp_ele() 
F = obj.calc_F(pd.R, pd.M)
F_λ_opt_ = obj.get_F_λ_opt(pd.M, pd.V, pd.S, pd.K, λ)

w = F_λ_opt_.optimizer 
F_λ_opt_.opt_val
F_λ_opt_.opt_stat
F_λ_opt_.λ
F_λ_opt_.sol_time


# pd_5 = stock_data.load_processed_data(1:5,1:5)
#-------------------------------------------- Projected gradient descent for finding optimizers-------------------------------------------- 
λ = spaces.gen_simp_ele() 
w     = obj.get_F_λ_optimizer_pgd(pd, λ, iter=10000)
w_acc = F_λ_optimal_pgd_acc = obj.get_F_λ_optimizer_pgd_acc(pd, λ, iter=10000)

F_λ_optimal_pgd     = obj.get_F_λ_optimal(w, λ, pd)
F_λ_optimal_pgd_acc = obj.get_F_λ_optimal(w_acc, λ, pd)
F_λ_optimal_IPOpt = obj.get_F_λ_optimal_IPOpt(pd, λ)


spaces.supp(w)
spaces.supp(w_acc)



hps =  spaces.get_λ_spaces(10, pd.R_max_max)
# fieldnames(typeof(λ_space))


F_λ_optimal_pgd, F_λ_optimal_IPOpt
abs((F_λ_optimal_pgd - F_λ_optimal_IPOpt)/F_λ_optimal_IPOpt)



λ_set = [spaces.gen_simp_ele(4,9)... ,rand(4)]
W = obj.get_F_λ_optimizer_pgd(pd, λ_set, iter=10000)
[spaces.supp(w) for w in W]

#-----------------------------------------------------------------------



w = obj.get_F_λ_optimizer_pgd(pd, λ,iter=10000)
spaces.is_simp(w)

∇F_λ = obj.get_∇F_λ(λ, pd.R, pd.M)
∇F_λ(w)

F_optimal_pgd   = obj.get_F_λ_optimal(w, λ, pd)
F_λ_optimal_IPOpt = obj.get_F_λ_optimal_IPOpt(pd, λ)
abs((F_optimal_pgd - F_λ_optimal_IPOpt)/F_λ_optimal_IPOpt)

#----------------------------------------------------------------------

# Testing IPOpt agianst gradient descent, unconstrained optimization of a convex function

## Test minₓ <x,r>² for some x
using JuMP
using Ipopt

n = 5
r = rand(n)
Q = r*r'

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "max_cpu_time", 60.0)
@variable(model, w[1:n])
@NLobjective(model, Min, sum( Q[i,j]*w[i]*w[j] for i in 1:n for j in 1:n))
optimize!(model)


## Grad. Desc.
include("pgd.jl")                   ; using .pgd

∇f =  w -> Q*w
x₀ = rand(n) / 10
Π = x -> x

w = proj_grad_dec(∇f, x₀, Π, 0.2, 10000) 
w_a = pgd.proj_acc_grad_dec(∇f, x₀, Π, 0.2, 10000) 

w'*Q*w, w_a'*Q*w_a, objective_value(model)



## Test minₓ <x,r>³ for some x
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "max_cpu_time", 60.0)
@variable(model, w[1:n])
@NLobjective(model, Min, sum( r[i]*r[j]*r[k]*w[i]*w[j]*w[k]  for i in 1:n for j in 1:n for k in 1:n))
optimize!(model)

∇f =  w -> 3*r*(r'*w)^2
w = pgd.proj_grad_dec(∇f, x₀, Π, 0.3, 10000)
w_acc = pgd.proj_acc_grad_dec(∇f, w, Π, 0.3, 10000)
∇f(w)

(r'*w)^3, objective_value(model)


x=y_prev=x₀
γ_seq =  pgd.γ(μ(10))
# for i ∈ 1:K
i = 2
β = 0.3
y = x - β*∇f(x)
x = (1-γ_seq[i])*y + γ_seq[i]*y_prev
y_prev = Π(y)
x = Π(x)
# end







