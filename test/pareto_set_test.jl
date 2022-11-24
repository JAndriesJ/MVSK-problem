module pareto_set_test
using Test    

scr_dir = dirname(pwd())*"/src/"
include(scr_dir*"stock_data.jl")     ; using .stock_data
include(scr_dir*"spaces.jl")         ; using .spaces
include(scr_dir*"pareto_set.jl")     ; using .pareto_set


pd  = stock_data.read_proc_data()
hps = spaces.get_λ_spaces(3, pd.R_max_max)

function run_tests()
    @testset "populate_results_csv" begin
        pareto_set.populate_results_csv(pd.R, pd.M, pd.V, pd.S, pd.K, hps.simp_λ_set, hps.simp_conv_mask; save_path="test_save.csv", box_size=0, sub=0)
        @test isfile("test_save.csv")
    end    

    @testset "read_results_csv" begin
        df = pareto_set.read_results_csv("test_save.csv")
        length(df.F1) == length(df.F3) == 10
        @test all(df.F2 .> 0)
        @test all(df.F4 .> 0)
        @test all(df.t .< 5)
        @test all(df.stat .== "LOCALLY_SOLVED")
        @test df.w[5][2] ≈ 0.04221480295912091
        @test df.is_conv == vec([1  1  1  0  0  0  1  1  0  1])
        @test sum.(df.λ) == ones(10) 

        rm("test_save.csv")
    end
end

end