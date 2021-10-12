module constraints
    include("stock_data.jl")
    include("moments.jl")
    using .moments ; const mom = moments
    using .stock_data

    export  make_PSD_constraint,
            make_localzing_constraints,
            make_localizing_ideal_constraint,
            make_objective_function

    """L([x]≦ₜ[x]ᵀ≦ₜ) ⪰ 0"""
    function make_PSD_constraint(N,t,Lx)
        Mₜ = mom.make_mon_expo(N,(t,t))
        return mom.get_Lxᵅ(Lx, Mₜ)
    end
    
    """L(xᵢ[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0"""
    function make_localzing_constraints(N,t,Lx)
        Mₜ₋₁ = mom.make_mon_expo(N,(t-1,t-1))
        return [mom.get_Lxᵅ(Lx, map(x -> x + mom.eᵢ(N,i), Mₜ₋₁)) for i ∈ 1:N]
    end
    
    """L((1 - ∑ᴺᵢ₌₁xᵢ)[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0"""
    function make_localizing_ideal_constraint(N,t,Lx)
        Mₜ₋₁ = mom.make_mon_expo(N,(t-1,t-1))
        return sum( [mom.get_Lxᵅ(Lx, map(x -> x + mom.eᵢ(N,i), Mₜ₋₁))  for i ∈ 1:N] ) 
    end

    """L(f(x)) = (-1)⁻¹*L(xᵀΦ⁽ᵏ⁾(x ⊗ ⋯ ⊗ x))"""
    function make_objective_function(N,k,Lx)
        dkm  = moments.make_mon_expo(N,k,isle=false)
        Lxk  = moments.get_Lxᵅ(Lx,dkm) 
        R    = load_relative_returns_data_matrix([1:N...])
        datk = map(α -> stock_data.get_data_moment(R,α),dkm)
        return (-1)^(k+1) * sum(Lxk .* datk)
    end
    
end