module grid_search
using Combinatorics, Test
export  Δ,
        min, 
        max

## Build a grid on the simplex
"""| {(x₁,...,xₖ) ∈ ℕᵏ  ~:~∑ᵢ xᵢ = n  }| = { n+k-1 \\choose k-1}""" 
function Δ_int(n,r) 
    r = r -1 # fix the logic here
    solutions = []
    for combo in [Combinatorics.combinations(0:(r+n-1),n-1)...]
        s = [combo[1]]
        for i in 2:(n-1)
            push!(s,combo[i]-combo[i-1]-1)
        end
        push!(s,r + n - 1 - combo[n - 1])
        push!(solutions,s)
    end
    return solutions
end
Δ(n,r) = Δ_int(n,r)./r

""" min f(x) s.t. x ∈ Δ(n,r) """
min(f,n,r) = minimum(f.(Δ(n,r)))

""" max f(x) s.t. x ∈ Δ(n,r) """
max(f,n,r) = maximum(f.(Δ(n,r)))

# Tests
function run_tests()
    @testset "Cardinality of set Δint_nr" begin
        for n ∈ rand(4:20,6), r ∈ rand(2:8,4)
            Δint_nr = Δ_int(n,r)
            @test size(unique(sum.(Δint_nr))) == (1,)
            @test unique(sum.(Δint_nr))[1] == r
            @test length(Δint_nr) == binomial(n+r-1, n-1)
        end
    end

    # @testset "max/min over grid" begin
    #     for n ∈ rand(4:20,6), r ∈ rand(2:8,4)
    #         f(x) = (rand(n,1)'*x)^rand(1:4,1)

    #         @test max(f,n,r) ==  max(f.(Δnr))
    #         @test min(f,n,r) ==  min(f.(Δnr))
    #     end
    # end

end  

end

