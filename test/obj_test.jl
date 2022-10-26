module obj_test
using Test    
using Random

include(pwd()*"/src/obj.jl")            ; using .obj
include(pwd()*"/src/spaces.jl")         ; using .spaces

Random.seed!(1234)

n = 10
R = reshape([1:n^2...],n,n)
M = [1:n...]

λ = spaces.gen_simp_ele()

function run_tests()
    @testset "Functions" begin
        F    = obj.calc_F(R,M)
        F_λ  = obj.calc_F_λ(λ,R,M)

        w = ones(n)
        @test [F[1](w), F[2](w), F[3](w), F[4](w)]  == [55, 255850.0,  1.300375e8, 6.6301333e10] 
        @test F_λ(w) == [-F[1](w), F[2](w), -F[3](w), F[4](w)]'*λ

    end    

    @testset "Optimizers" begin

    end
end

end


## Tests
# w = rand(20)
# f₁(w), f₂(w) ,f₃(w) ,f₄(w)
# ∇f₁(), ∇f₂(w), ∇f₃(w), ∇f₄(w)
# F(λ,w)
# ∇F(w) = ∇F(λ,w)

# F_mult(w) = [f₁(w), f₂(w), f₃(w), f₄(w)]

# function get_F_λ_optimizers(λ_set, pd; simp_check=false, Π=Π, β=β, iter=iter)
#     if simp_check
#         return map(λ -> get_F_λ_optimizer(λ, pd.R, pd.M, Π=Π, β=β, iter=iter), λ_set)
#     else
#         return map( λ -> spaces.is_simp(λ) ? get_F_λ_optimizer(λ, R, M,  Π=Π, β=β, iter=iter) : NaN, λ_set)
#     end
# end

# w_star = obj.get_F_λ_optimizer(λ, R, M, Π, β, iter)
        # @test w_star == [1.0, zeros(9)...]

        # λ_set = [spaces.Euc_proj_simp(rand(4)) for i in 1:10]
        # w_star_set = obj.get_F_λ_optimizers(λ_set, R, M, false, Π, β, iter)
        # @test w_star_set == [9.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0]

        # # obj.get_F_optimals(w_star_set, R, M, tru