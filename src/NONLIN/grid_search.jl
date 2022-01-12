module grid_search
using Combinatorics,  DataFrames,Test
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

""" min f(x) s.t. x ∈ Δ(n,r) , argmin """
function min_Δ(f,n,r)
    Δnr = Δ(n,r)
    fΔnr = f.(Δnr)
    min_f = minimum(fΔnr)
    return Δnr[fΔnr .== min_f][1], min_f
end


""" max f(x) s.t. x ∈ Δ(n,r) """
function max_Δ(f,n,r)
    Δnr = Δ(n,r)
    fΔnr = f.(Δnr)
    max_f = maximum(fΔnr)
    return Δnr[fΔnr .== max_f][1], max_f
end

function batch_min(f,n,r_range)
    df = ini_df()    
    for r in r_range
        println("Now running gridsearch with density: $r")
        sol_t = @elapsed  obj_arg, obj_val = min_Δ(f,n,r)
        push!(df,  (r, obj_val, obj_arg, sol_t))
    end
    return df
end

function batch_max(f,n,r_range)
    df = ini_df()         
    for r in r_range
        println("Now running gridsearch with density: $r")
        sol_t = @elapsed  obj_arg, obj_val = max_Δ(f,n,r)
        push!(df,  (r, obj_val, obj_arg, round(sol_t)))
    end
    return df
end

function ini_df()
    return  DataFrame(  grid_density=Int[], 
                        obj_val=Float64[],
                        obj_arg=Vector{Float64}[],
                        time_s=Float64[])  
end


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

    # val_1, f_val_1 = grid_search.max_Δ(f₁, mod.n_stocks,30)
    # val_2, f_val_2 = grid_search.min_Δ(f₂, mod.n_stocks,30)
    # val_3, f_val_3 = grid_search.max_Δ(f₃, mod.n_stocks,30)
    # val_4, f_val_4 = grid_search.min_Δ(f₄, mod.n_stocks,30)
    # f₂(val_2) ==  f_val_2


end  

end

