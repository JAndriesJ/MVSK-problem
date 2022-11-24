module spaces_test
using Test
using Random

scr_dir = dirname(pwd())*"/src/"
include(scr_dir*"spaces.jl") ; using .spaces

Random.seed!(1234)
function run_tests()
    @testset "simplex" begin
        Π = x -> spaces.Euc_proj_simp(x)
        W = map(x -> Π(x), [rand(5) for i ∈ 1:10])
        @test all([spaces.is_simp(w) for w ∈ W])
        @test all([sum(w) ≈ 1.0 for w ∈ W])
    end

    @testset "gen_simp_ele" begin
        w = spaces.gen_simp_ele()
        @test length(w) == 4
        @test sum(w) ≈ 1
        @test all(w .≥ 0) 

        n = rand(5:10)
        w = spaces.gen_simp_ele(n)   
        @test length(w) == n  
        @test sum(w) ≈ 1
        @test all(w .≥ 0) 
        
        k = rand(10:30)
        W = spaces.gen_simp_ele(n,k) 
        length(W) == k
        for w in W
            @test length(w) == n  
            @test sum(w) ≈ 1
            @test all(w .≥ 0)
        end
    end 

    @testset "r_simplex" begin    
        W = map(x -> spaces.Euc_proj_simp(x)[1:4], [rand(5) for i ∈ 1:10])
        @test all([spaces.is_r_simp(w) for w ∈ W])
    end

    @testset "sphere" begin
        N = x -> spaces.L2_normalize(x)
        W = map(x -> N(x), [rand(5) for i ∈ 1:10])
        @test all([spaces.is_sphere(w) for w ∈ W])
        @test all([sum(w.^2) ≈ 1.0 for w ∈ W])
    end

    @testset  "get_λ_spaces" begin
        mf = 10
        b = 0.4433
        λ_spaces = spaces.get_λ_spaces(mf, b)
        @test fieldnames(typeof(λ_spaces)) == (:tics, :λ2, :λ3, :λ4, :λ, :simp_λ_set, :simp_mask, :conv_mask, :simp_conv_mask, :conv_rel_vol)
        @test all(sum.(λ_spaces.simp_λ_set) .≈ ones(sum(λ_spaces.simp_mask)))
        @test sum(λ_spaces.simp_mask .+ λ_spaces.conv_mask .== 2) == sum(λ_spaces.simp_conv_mask)
    end

end

end
