module Moments
using Test
export  eᵢ,
        make_mon_expo,
        get_Lxᵅ,
        expo_kron,
        run_tests

"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [Int(j==i) for j in 1:n]

""" [x]≦ₜ or [x]₌ₜ """
function make_mon_expo(n::Int,t::Int; isle::Bool = true)
    t == 0 ? (return [eᵢ(n,0)]) : 0
    mon_expo = make_mon_expo(n,t-1;isle=isle)
    mon_expo_vec = reshape([m + eᵢ(n,i) for i ∈ 1:n, m ∈ mon_expo],:,1)
    return unique(isle ? vcat(mon_expo, mon_expo_vec) : mon_expo_vec)
end

""" [x]≦ₜ([x]≦ₜ)ᵀ or [x]₌ₜ([x]₌ₜ)ᵀ """
function make_mon_expo(n::Int,t::Tuple{Int,Int}; isle::Bool = true)
    mon_expo_vec1      = make_mon_expo(n,t[1]; isle=isle)
    mon_expo_vec2      = make_mon_expo(n,t[2]; isle=isle)
    return [mi+mj for mi in mon_expo_vec1, mj in mon_expo_vec2]
end

""" (Lx,α) -> var[x^α] """
get_Lxᵅ(Lx,α) = isempty(α) ?  0.0 : map(y -> Lx[y],α)

"""
A = [αᵢⱼ]ᵢ,ⱼ ∈ (ℕⁿ)ᴺˣᴹ ; B = [βₖₗ]ₖ,ₗ ∈ (ℕⁿ)ᴷˣᴸ
expo_kron(A,B) = [αᵢⱼ+βₖₗ]ᵢₖ,ⱼₗ ∈ (ℕⁿ)ᴺᴷˣᴹᴸ
"""
function expo_kron(A,B)
    C_temp = [ a + b for  b in B, a ∈ A]
    rcs =  [(i,j) for i in 1:size(A)[1] for j in 1:size(B)[1]]
    return [C_temp[ij[2],kl[2],ij[1],kl[1]] for ij in rcs , kl in rcs ]
 end

function run_tests()
    ## moments

    @testset "eᵢ" begin
        n = 7
        for k in 1:n
            e = eᵢ(n,k)
            @test length(e) == n
            @test maximum(e) == 1
            @test minimum(e) == 0
            @test sum(e) == 1
        end
    end

    @testset "make_mom_expo_mat" begin
        mon_expo_mat = make_mon_expo(3,(1,1);isle = true)
        @test [mon_expo_mat[i,i] for i in 1:4]  == [ [0, 0, 0],[2, 0, 0],[0, 2, 0],[0, 0, 2]]
        n,t = 4,2
        M_vec =  make_mon_expo(n,t)
        @test length(M_vec) == binomial(n+t,t)
        M_mat =  make_mon_expo(n,(t,t))
        @test M_vec == M_mat[1,:]
    end

    # @testset "(Lx,α) -> var[x^α]" begin
    #     N,t = 5,4
    #     model = Model()
    #     moment_vec = make_mon_expo(N,2*t)
    #     @variable(model, Lx[moment_vec]) # Create variables
    #     α_list  = [[0, 1, 0, 0, 0], [0, 0, 0, 4, 4], [0, 0, 0, 0, 8]]
    #     @test get_Lxᵅ(Lx,α_list) == [Lx[[0, 1, 0, 0, 0]], Lx[[0, 0, 0, 4, 4]], Lx[[0, 0, 0, 0, 8]]]
    #     β_list  =  [[0, 6, 0, 1, 0],[3, 3, 0, 0, 2],[0, 2, 6, 0, 0],[2, 0, 0, 1, 3]]
    #     @test get_Lxᵅ(Lx,α_list) == [Lx[[0, 6, 0, 1, 0]], Lx[[3, 3, 0, 0, 2]], Lx[[0, 2, 6, 0, 0]], Lx[[2, 0, 0, 1, 3]]]
    # end
end


end




## Cardinality restricted
# """[̃x]≦ₜ or [̃x]₌ₜ"""
# function get_cardinality_restricted_mon_expo(n::Int,t::Int,ℓ::Int; isle::Bool = true)
#     mon_expo = make_mon_expo(n,t; isle)
#     return mon_expo[[sum(.!iszero.(expo)) < ℓ  for expo in mon_expo]]
# end

# """ [̃x]≦ₜ([̃x]≦ₜ)ᵀ or [̃x]₌ₜ([̃x]₌ₜ)ᵀ """
# function get_cardinality_restricted_mon_expo(n::Int,t::Tuple{Int,Int},ℓ::Int; isle::Bool = true)
#     mon_expo_vec1      = get_cardinality_restricted_monomials(n,t[1],ℓ; isle=isle)
#     mon_expo_vec2      = get_cardinality_restricted_monomials(n,t[2],ℓ; isle=isle)
#     return [mi+mj for mi in mon_expo_vec1, mj in mon_expo_vec2]
# end