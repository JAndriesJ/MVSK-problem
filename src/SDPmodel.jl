module SDPmodel
    using JuMP
    include("moments.jl")
    include("SDPconstraints.jl")
    using .moments 
    using .SDPconstraints 

    export get_SDP_model

    function get_SDP_model(N,t,k)
        model = Model()
        @variable(model, Lx[moments.make_mon_expo(N,2*t)]) # Create variables
        ## L(1) = 1
        @constraint(model, Lx[moments.make_mon_expo(N,0)[1,1]] == 1)
        ## L([x]≦ₜ[x]ᵀ≦ₜ) ⪰ 0
        @constraint(model, SDPconstraints.make_PSD_constraint(N,t,Lx) in PSDCone())
        ## L(xᵢ[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0 ∀ i ∈ [N]
        for con in SDPconstraints.make_localzing_constraints(N,t,Lx)
            @constraint(model, con in PSDCone())
        end
        ## L((1 - ∑ᴺᵢ₌₁xᵢ)[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0
        loc_ideal_con = SDPconstraints.make_localizing_ideal_constraint(N,t,Lx)
        @constraint(model, loc_ideal_con .== zeros(size(loc_ideal_con)))
        ## Min f
        f = SDPconstraints.make_objective_function(N,k,Lx)
        @objective(model, Min, f) #  Set objective
        return model 
    end

end