module spaces

using PlotlyJS
include("λ_char.jl")                ; using .λ_char 

export get_λ_spaces,
       is_simp,
       is_r_simp,
       gen_simp_ele,
       Euc_proj_simp,
       is_sphere,
       L2_normalize,
       supp
  
#-----------------------------------------------------Simplex-----------------------------------------------------
struct hyper_parameter_space
    tics
    λ2
    λ3
    λ4
    λ
    simp_λ_set
    simp_mask
    conv_mask 
    simp_conv_mask 
    conv_rel_vol
end

function get_λ_spaces(mf, b)
    r                    = range(0, stop=1, length=mf)
    λ2, λ3, λ4           = mgrid(r, r, r)
    λ_grid               = gen_λ_grid(λ2,λ3,λ4,mf)
    simp_mask, conv_mask = calc_simp_masks(λ_grid, b)
    
    simp_λ_set     = λ_grid[simp_mask]
    simp_conv_mask = conv_mask[simp_mask] 
    conv_rel_vol   = sum(simp_mask[conv_mask])

    return hyper_parameter_space(   mf,
                                    λ2,
                                    λ3,
                                    λ4,
                                    λ_grid,
                                    simp_λ_set,
                                    simp_mask,
                                    conv_mask,
                                    simp_conv_mask,
                                    conv_rel_vol)
end
gen_λ_grid(λ2,λ3,λ4,mf) = reshape([[1-(λ2[i]+λ3[i]+λ4[i]), λ2[i], λ3[i], λ4[i]] for i ∈ 1:mf^3],(mf, mf, mf))
function calc_simp_masks(λ_grid, b)
    temp = [[is_simp(λ) ? true : false , λ_char.λ₃_bound_Δ(λ, b) ? true : false]   for λ ∈ λ_grid]

    simp_mask = map(x -> x[1], temp)
    conv_mask = map(x -> x[2], temp)
    return  simp_mask, conv_mask
end

is_simp(λ) = sum(λ) ≈ 1 && all(λ .≥ 0)
is_r_simp(λ) =  sum(λ) ≤ 1  && all(λ .≥ 0)

gen_simp_ele(n)   =  Euc_proj_simp(rand(n))
gen_simp_ele()    =  gen_simp_ele(4)
gen_simp_ele(n,k) = [gen_simp_ele(n)  for  i in 1:k] 

function Euc_proj_simp(x)
    u = sort(x,rev=true)
    p = maximum(calc_val_vec(u))
    λ = 1/p*(1-sum(u[1:p]))
    return [maximum([x[i] + λ, 0]) for i ∈ 1:length(u)]
end
calc_val_vec(u) = [j for j ∈ 1:length(u) if u[j] + 1/j*(1-sum(u[1:j])) > 0]

#-----------------------------------------------------Sphere-----------------------------------------------------   
L2_norm(x) = sqrt(sum(x.^2))
L2_normalize(x) = x/L2_norm(x)
is_sphere(x) =L2_norm(x) ≈  1
#----------------------------------------------------Box-----------------------------------------------------   



#-----------------------------------------------------Utility-----------------------------------------------------   
supp(w) = sum(w .> 0)

end


# function ext_λ(hps)
#     sm = hps.simp_mask 
#     simp_λ_set = hps.λ[sm]
#     simp_conv_mask = hps.conv_mask[sm] 
#     return simp_λ_set, simp_conv_mask 
# end