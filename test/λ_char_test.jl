module λ_char_test
using Test
using Random

include(pwd()*"/src/spaces.jl")         ; using .spaces
include(pwd()*"/src/λ_char.jl")         ; using .λ_char

Random.seed!(1234)

function run_tests()
    @testset "λ₃_bounds" begin
        r = 1.5
        Λ = spaces.gen_simp_ele(3,10)
        @test [λ_char.λ₃_bound(λ[1],λ[2],λ[2]) for λ ∈ Λ] == [1,  0 , 1,  1,  0,  1,  0,  1,  1,  1]
        @test [λ_char.λ₃_bound_Δ(λ[1],λ[2],λ[2],r) for λ ∈ Λ] == [1,  1 , 1,  1,  1,  1,  1,  1,  1,  1]
       
       for _ in 1:7
            λ = spaces.gen_simp_ele(3)
            r = rand(1)[1] + 1.0
            @test λ_char.λ₃_bound(λ[1], λ[2]) ≤ λ_char.λ₃_bound_Δ(λ[1], λ[2], r)
       end
    end 
end

end