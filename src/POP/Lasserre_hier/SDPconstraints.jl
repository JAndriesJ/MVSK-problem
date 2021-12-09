module SDPconstraints
    include("moments.jl")
    using .Moments ; const mom = Moments

    export  make_variables,
            make_L1,
            make_PSD_constraint,
            make_localzing_constraints,
            make_localizing_ideal_constraint,
            make_objective_function

    # TODO: Make N and Lx global variable in this moduel

    """L([x]≦₂ₜ)"""
    make_variables(N,t) = mom.make_mon_expo(N,2*t)

    """L(1)"""
    make_L1(N,Lx) = Lx[mom.make_mon_expo(N,0)[1,1]]

    """L([x]≦ₜ[x]ᵀ≦ₜ) ⪰ 0"""
    function make_PSD_constraint(N,Lx,t)
        Mₜ = mom.make_mon_expo(N,(t,t))
        return mom.get_Lxᵅ(Lx, Mₜ)
    end

    """L(xᵢ[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0""" # TODO: MAKE THIS STRONGER BY INCLUDING THE DEGREE 2t terms
    function make_localzing_constraints(N,Lx,t)
        Mₜ₋₁ = mom.make_mon_expo(N,(t-1,t-1))
        return [mom.get_Lxᵅ(Lx, map(x -> x + mom.eᵢ(N,i), Mₜ₋₁)) for i ∈ 1:N]
    end
    
    """L((1 - ∑ᴺᵢ₌₁xᵢ)[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0""" # TODO: MAKE THIS STRONGER BY INCLUDING THE DEGREE 2t terms
    function make_localizing_ideal_constraint(N,Lx,t)
        Mₜ₋₁ = mom.make_mon_expo(N,(t-1,t-1))
        return make_PSD_constraint(N,Lx,t-1) - sum( [mom.get_Lxᵅ(Lx, map(x -> x + mom.eᵢ(N,i), Mₜ₋₁))  for i ∈ 1:N] ) 
    end

    """L(f(x)) = L(xᵀ skewness_matrix (x ⊗ x)) or L((x ⊗ x)ᵀ kurtosiss_matrix (x ⊗ x))"""
    function make_objective_function(N,Lx,Data_matrix) # Flawed.
        n,m = size(Data_matrix); n==m ? k = 4 : k = 3
        MomentMat_expo_deg_eq_k  = mom.make_mon_expo(N,k,isle=false)
        MomentMat_deg_eq_k  = mom.get_Lxᵅ(Lx,MomentMat_expo_deg_eq_k) 
        println(size(MomentMat_deg_eq_k))
        return sum(MomentMat_deg_eq_k .* Data_matrix) ## f(x) = ∑ L(xᵅ) Data_matrixᵅ
    end
    
end

  # R    = stock_data.load_relative_returns_data_matrix([1:N...])
    # datk = map(α -> stock_data.get_data_moment(R,α),dkm) ## ???
    
    # function run_tests()
    #     using Test, JuMP
    #     ## constraints
    #     @testset "all the constraints" begin
    #         model = JuMP.Model()
    #         N,t = 5,3
    #         @variable(model, Lx[moments.make_mon_expo(N,2*t)]) # Create variables

    #         @test Lx[moments.make_mon_expo(N,0)[1,1]] == Lx[[0, 0, 0, 0, 0]]   ## L(1) 
    #         PSD_moment_matrix     = constraints.make_PSD_constraint(N,t,Lx) ## L([x]≦ₜ[x]ᵀ≦ₜ)
    #         @test size(PSD_moment_matrix) == (binomial(N+t,t),binomial(N+t,t))
    #         @test PSD_moment_matrix[end,end] == Lx[moments.eᵢ(N,N)*2*t]

    #         localzing_constraints = [con for con in constraints.make_localzing_constraints(N,t,Lx)] ## L(xᵢ[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0 ∀ i ∈ [N]
    #         @test length(localzing_constraints) == N
    #         @test size(localzing_constraints[1]) == (binomial(N+t-1,t-1),binomial(N+t-1,t-1))
    #         @test localzing_constraints[2][end,end] == Lx[[0, 1, 0, 0, 4]]

    #         loc_ideal_con         = constraints.make_localizing_ideal_constraint(N,t,Lx) ## L((1 - ∑ᴺᵢ₌₁xᵢ)[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0
    #         @test size(loc_ideal_con) == (binomial(N+t-1,t-1),binomial(N+t-1,t-1))
    #         @test loc_ideal_con[end,end] == Lx[[0, 0, 0, 0, 4]] - Lx[[1, 0, 0, 0, 4]] - Lx[[0, 1, 0, 0, 4]] - Lx[[0, 0, 1, 0, 4]] - Lx[[0, 0, 0, 1, 4]] - Lx[[0, 0, 0, 0, 5]]

    #         φ⁽³⁾ = constraints.make_objective_function(N,3,Lx) ## wᵀΦ⁽³⁾(w⊗w)
    #         φ⁽⁴⁾ = constraints.make_objective_function(N,4,Lx) ## wᵀΦ⁽⁴⁾(w⊗w⊗w)
    #         @test φ⁽³⁾.terms.keys == moments.get_Lxᵅ(Lx,moments.make_mon_expo(N,3,isle=false) )
    #         @test φ⁽⁴⁾.terms.keys == moments.get_Lxᵅ(Lx,moments.make_mon_expo(N,4,isle=false) )
    #     end
    # end