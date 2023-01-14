module pareto_set

using DataFrames, CSV
using Dates
using LinearAlgebra

include("obj.jl")                   ; using .obj

export  prep_results_csv,
        read_results_csv
 
function create_empty_results_csv(save_path = pwd()*"\\"*replace(string(now()),":"=>"-")*s*".csv")
    if isfile(save_path)
        println("Using provided save path $save_path")
    else
        touch(save_path)
        println("Creating save path $save_path")
        open(save_path, "w") do io
            write(io, "λ|is_conv|w|F1|F2|F3|F4|Fλ|t|stat\n")
        end;
    end
    
    return save_path
end

function prep_results_csv(pd,λ_set,conv_mask, save_path::String; iter::Int=0 ,box_size=0, sub=0)
    isfile(save_path) ? nothing : create_empty_results_csv(save_path)
    F = obj.calc_F(pd.M, pd.V, pd.S_mat, pd.K_mat)
    populate_results_csv(F, pd, λ_set, conv_mask, save_path, iter, box_size, sub)
    return save_path
end

function populate_results_csv(F, pd, λ_set, conv_mask, save_path, iter, box_size, sub)
    X = []
    n = length(pd.M)
    np = length(λ_set)
    for i ∈ 1:np 
        λ = λ_set[i]    
        x₀ =  get_near_x₀(X,n,λ_set,λ,i) 

        F_λ_opt_ =  obj.get_F_λ_opt(pd,λ,x₀=x₀,iter=iter, box_size=box_size, sub=sub)
        w, st, os, Fw, F_λw =  proc_results_line(F_λ_opt_, F,λ)

        push!(X,w)
        write_results(save_path, λ, conv_mask[i], w, Fw, F_λw, st, os)
    end
end
get_near_x₀(X,n,λ_set,λ,i) =  isempty(X) ? zeros(n) : X[argmin(map(x -> norm(x-λ), λ_set[1:i-1]))]


function proc_results_line(F_λ_opt_,F,λ)
    w        = F_λ_opt_.optimizer 
    st       = F_λ_opt_.sol_time
    os       = F_λ_opt_.opt_stat
    Fw       =  F[1](w),F[2](w),F[3](w),F[4](w)
    F_λw     = [-λ[1],λ[2],-λ[3],λ[4]]'*[F[1](w),F[2](w),F[3](w),F[4](w)]
    return w, st, os, Fw, F_λw
end

read_results_csv(csv_path, delim="|") = clean_results_csv(
                                               CSV.read(csv_path, DataFrame, delim=delim)
                                                )
function clean_results_csv(df)
    df.λ       = txt_to_math.(df.λ) 
    df.w       = txt_to_math.(df.w)
    return df
end
txt_to_math(s) = eval(Meta.parse(s))

function write_results(save_path,λ,is_con,w,Fw,F_λw,st,os)
    open(save_path, "a") do io
        write(io, "$(λ)|$(is_con)|$(w)|$(Fw[1])|$(Fw[2])|$(Fw[3])|$(Fw[4])|$(F_λw)|$(st)|$(os)\n")
    end
end


end







