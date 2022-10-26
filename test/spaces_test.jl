module spaces_test
using Test

include(pwd()*"/src/spaces.jl")         ; using .spaces

function run_tests()
    @testset "simplex" begin
        Π = x -> spaces.Euc_proj_simp(x)
        W = map(x -> Π(x), [rand(5) for i ∈ 1:10])
        @test all([spaces.is_simp(w) for w ∈ W])
        @test all([sum(w) ≈ 1.0 for w ∈ W])
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
end

end


# get_λ_spaces,
# is_simp,
# is_r_simp,
# gen_simp_ele,
# Euc_proj_simp,
# is_sphere,
# L2_normalize,
# supp