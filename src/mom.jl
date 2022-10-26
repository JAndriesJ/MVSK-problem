module mom
export  eᵢ,
        make_mon_expo

"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [Int(j==i) for j in 1:n]
#make_mon_expo-------------------------------------------------------------------------------------------------------------------------------------
"""[x]≦ₜ := [xᵅ for all α ∈ ℕⁿₜ] or [x]₌ₜ  := [xᵅ for all α ∈ ℕⁿ≤ₜ]"""
function make_mon_expo(n::Int,t::Int; isle::Bool = true)
    t == 0 ? (return [eᵢ(n,0)]) : 0
    tmp = make_mon_expo(n,t-1;isle=isle)
    M_vec = reshape([m + eᵢ(n,i) for i ∈ 1:n, m ∈ tmp],:,1)
    return unique(isle ? vcat(tmp, M_vec) : M_vec)
end

end
