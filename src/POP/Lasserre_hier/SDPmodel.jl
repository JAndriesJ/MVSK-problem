module SDPmodel
    using JuMP
    include("SDPconstraints.jl")
    using .SDPconstraints 

    export get_SDP_model

    function get_SDP_model(Data_matrix,t)
        n,m = size(Data_matrix)
        @assert m == n^2 || m == n
        n == m ? N = Int(sqrt(n)) : N = n
        
        model = Model()
        @variable(model, Lx[SDPconstraints.make_variables(N,t)]) 
        ## L(1) = 1
        @constraint(model, SDPconstraints.make_L1(N,Lx) == 1)
        ## L([x]≦ₜ[x]ᵀ≦ₜ) ⪰ 0
        @constraint(model, SDPconstraints.make_PSD_constraint(N,Lx,t) in PSDCone())
        ## L(xᵢ[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0 ∀ i ∈ [N]
        for con in SDPconstraints.make_localzing_constraints(N,Lx,t)
            @constraint(model, con in PSDCone())
        end
        ## L((1 - ∑ᴺᵢ₌₁xᵢ)[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0
        loc_ideal_con = SDPconstraints.make_localizing_ideal_constraint(N,Lx,t)
        @constraint(model, loc_ideal_con .== zeros(size(loc_ideal_con)))
        ## Min f
        f = SDPconstraints.make_objective_function(N,Lx,Data_matrix)
        @objective(model, Min, f) 
        return model 
    end

end