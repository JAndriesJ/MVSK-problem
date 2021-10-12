module SDPoptimized
export optimize_SDP

using MosekTools, JuMP
# Pkg.build("MosekTools")  # [[[In the powershell]]] # Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy Unrestricted 
function optimize_SDP(model)
    JuMP.set_optimizer(model, Mosek.Optimizer)
    # optimize
    optimize!(model)

    # output results
    println("Primal: ", primal_status(model))
    println("Dual: ", dual_status(model))
    println("Objective: ", objective_value(model))
    return model
end
    
end