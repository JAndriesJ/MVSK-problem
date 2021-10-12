module SDPoptimized
include("src\\SDPmodel.jl")
using .SDPmodel 


export optimize_SDP,
       batch_optimize_SDP

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

using Dates

function batch_optimize_SDP(N_list,t_list,k)
    save_dir = "C:\\Users\\jandr\\MahCodes\\MVSK\\assets\\$(string(Dates.today()))-bounds_data.md"
    abc = open(save_dir, "w")
    write(abc, "## φ^($(k)) \n \n" )
    write(abc,"|level|#stocks|primal_status|dual_status|obj_val|comp_time(s)|\n")
    write(abc,"|---|---|---|---|\n")
    close(abc)

    for t in t_list
        for N in N_list 
            # try 
            SDPmod = SDPmodel.get_SDP_model(N,t,4)
            optSDPmod = SDPoptimized.optimize_SDP(SDPmod)
            P_stat = string(primal_status(optSDPmod)) 
            D_stat = string(dual_status(optSDPmod))
            obj_val = JuMP.objective_value(optSDPmod)
            comp_time = JuMP.solve_time(optSDPmod)
            abc = open(save_dir, "a")
                write(abc,"|$(t)|$(N)|$(P_stat)|$(D_stat)|$(obj_val)|$(comp_time)|\n")
            close(abc)
            # catch 
            # abc = open(save_dir, "a")
            #     write(abc,"|$(t)|$(N)|Mem_err.|-|-|-|\n")
            # close(abc)
            # break
            # end
        end
    end
end


    
end