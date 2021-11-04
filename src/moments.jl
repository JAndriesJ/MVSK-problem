module moments

export  eᵢ,
        make_mon_expo,
        get_Lxᵅ

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

## Cardinality restricted

"""[̃x]≦ₜ or [̃x]₌ₜ"""
function get_cardinality_restricted_mon_expo(n::Int,t::Int,ℓ::Int; isle::Bool = true)
    mon_expo = make_mon_expo(n,t; isle)
    return mon_expo[[sum(.!iszero.(expo)) < ℓ  for expo in mon_expo]]
end

""" [̃x]≦ₜ([̃x]≦ₜ)ᵀ or [̃x]₌ₜ([̃x]₌ₜ)ᵀ """
function get_cardinality_restricted_mon_expo(n::Int,t::Tuple{Int,Int},ℓ::Int; isle::Bool = true)
    mon_expo_vec1      = get_cardinality_restricted_monomials(n,t[1],ℓ; isle=isle)
    mon_expo_vec2      = get_cardinality_restricted_monomials(n,t[2],ℓ; isle=isle)
    return [mi+mj for mi in mon_expo_vec1, mj in mon_expo_vec2]
end

""" (Lx,α) -> var[x^α] """
get_Lxᵅ(Lx,α) = isempty(α) ?  0.0 : map(y -> Lx[y],α)

end
