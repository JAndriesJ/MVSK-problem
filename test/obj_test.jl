module obj_test
using Test    

main_dir = dirname(pwd())
println(main_dir)

include(main_dir*"\\src\\"*"stock_data.jl")   ; using .stock_data 
include(main_dir*"/src/"*"spaces.jl")         ; using .spaces
include(main_dir*"/src/"*"obj.jl")            ; using .obj

global pd = stock_data.get_proc_data(main_dir*"\\assets\\stock_prices.csv")

function run_tests()
    @testset "Functions" begin
        F    = obj.calc_F(pd.R, pd.M)
        w = ones(length(pd.M))
        @test F[1](w) ≈ sum(pd.M)
        @test F[2](w) ≈ w'*pd.V*w
        # @test F[3](w) == -0.5903949337088907
        # @test F[4](w) ==  0.5468435597841993
    end    

    @testset "Optimizers" begin
        λ = [0.2878  0.302  0.2831  0.1271]
        # simplex 
        F_λ_opt =  obj.get_F_λ_opt(pd.M, pd.V, pd.S, pd.K, λ, box_size=0, sub=0);
        @test F_λ_opt.domain == "simplex"
        #@test all(F_λ_opt.optimizer[1:2] .≈ [1.3754765623294924e-5, 1.1595362317398075e-5])
        #@test F_λ_opt.opt_val ≈ -0.00023835052763281075
        @test F_λ_opt.opt_stat == "LOCALLY_SOLVED"
        #@test isapprox(F_λ_opt.sol_time, 1.1430001258850098, atol =1)
        @test F_λ_opt.λ == λ 
        # unit box
        F_λ_opt =  obj.get_F_λ_opt(pd.M, pd.V, pd.S, pd.K, λ, box_size=1, sub=0);
        @test F_λ_opt.domain == "box"
        #@test all(F_λ_opt.optimizer[1:3]' .≈ [ -0.9993097208151496 -0.7954362484571216 -0.1983273758207766])
        #@test F_λ_opt.opt_val ≈ -0.0008935698828165571
        @test F_λ_opt.opt_stat == "LOCALLY_SOLVED"
        #@test isapprox(F_λ_opt.sol_time, 1.1430001258850098, atol =1)
        @test F_λ_opt.λ == λ  
        # unit sparse   

    end
end

end
