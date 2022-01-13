module MVSK
using CSV, DataFrames, Dates

include("mod\\stock_data.jl")   ; using .stock_data ; const sd = stock_data
include("mod\\data_moments.jl") ; using .data_moments ; const dm = data_moments
include("mod\\objective_functions.jl") ; using .objective_functions ; const obf = objective_functions

include("opt\\SOS.jl") ; using .SOS
include("opt\\grid_search.jl") ; using .grid_search  ; const gs = grid_search
include("opt\\convex_quadratic_optimization.jl") ; using .convex_quadratic_optimization ; const cqo = convex_quadratic_optimization  

export get_MVSK_data,
       get_MVSK_bounds,
       save_MVSK

mutable struct MVSK_data
    n_stocks::Int64
    Data_matrix::Matrix{Float64}
    moment_matrices::Dict{String, Matrix{Float64}}
    #stock_names::Vector{String}
    # begin_date::Dates.Date
    # end_date::Dates.Date
end

mutable struct MVSK_bounds
    λ::Vector{Float64}
    deg_range::StepRange{Int64, Int64}
    grid_range::Vector{Int64}
    f₁ₘₐₓ::Dict
    f₂ₘᵢₙ::Dict
    f₃ₘₐₓ::Dict
    f₄ₘᵢₙ::Dict
    𝒻_tilde::Dict
end    

###################################################################################
function get_MVSK_data(n_stocks = 6)
    Data_matrix = sd.load_cont_comp_returns()[:,1:n_stocks]
    moment_matrices = Dict( "means_vector"      => dm.get_means_vector(Data_matrix), 
                            "covariance_matrix" => dm.get_covariance_matrix(Data_matrix),
                            "skewness_matrix"   => dm.get_skewness_matrix(Data_matrix),
                            "kurtosis_matrix"   => dm.get_kurtosis_matrix(Data_matrix))
    
    MVSK_d = MVSK_data(n_stocks,
                       Data_matrix,
                       moment_matrices)

    return MVSK_d
end

function get_MVSK_bounds(mod, λ = ones(4)*0.25, deg_range = 4:2:6, grid_range = 2 .^ (1:5))
    mom_mats = mod.moment_matrices
    f₁ₘₐₓ = get_f₁ₘₐₓ(mom_mats["means_vector"])
    f₂ₘᵢₙ = get_f₂ₘᵢₙ(mom_mats["covariance_matrix"])
    f₃ₘₐₓ = get_f₃ₘₐₓ(mom_mats["skewness_matrix"], deg_range, grid_range)            
    f₄ₘᵢₙ = get_f₄ₘᵢₙ(mom_mats["kurtosis_matrix"], deg_range, grid_range)
    𝒻_tilde = get_𝒻_tilde(λ, mom_mats, f₁ₘₐₓ, f₂ₘᵢₙ, f₃ₘₐₓ, f₄ₘᵢₙ, grid_range)
  
    MVSK_b =    MVSK_bounds(λ,
                            deg_range,
                            grid_range,
                            f₁ₘₐₓ,f₂ₘᵢₙ,f₃ₘₐₓ,f₄ₘᵢₙ,
                            𝒻_tilde
                            )

    return MVSK_b
end

###################################################################################
function get_f₁ₘₐₓ(me_vec)
    f₁ₘₐₓ_val, f₁ₘₐₓ_arg = maximum(me_vec), me_vec .== maximum(me_vec)
    return Dict("f₁ₘₐₓ_val"=> f₁ₘₐₓ_val,
                "f₁ₘₐₓ_arg"=> f₁ₘₐₓ_arg)
end

function get_f₂ₘᵢₙ(cov_mat)
    f₂ₘᵢₙ_val = cqo.get_variance_min(cov_mat*(10^4))*(10^(-4)) # scaling for stability
    return Dict("f₂ₘᵢₙ_val"=> f₂ₘᵢₙ_val)
end

function get_f₃ₘₐₓ(sk_mat, deg_range, grid_range)
    n_stocks = size(sk_mat)[1] 
    f₃ = obf.get_f₃(sk_mat) 
    
    f₃ₘₐₓ_pop_df = SOS.get_SOS_skewness_estimates(sk_mat, deg_range)
    f₃ₘₐₓ_gs_df = gs.batch_max(f₃, n_stocks, grid_range)
    return Dict("SOS_estimates"=> f₃ₘₐₓ_pop_df,
                "grid_search_estimates"=> f₃ₘₐₓ_gs_df)
end

function get_f₄ₘᵢₙ(kur_mat, deg_range, grid_range)
    f₄ = obf.get_f₄(kur_mat)
    n_stocks = Int(sqrt(size(kur_mat)[1]))

    f₄ₘᵢₙ_cqr = cqo.get_relaxed_kurtosis(kur_mat) 
    f₄ₘᵢₙ_pop_df = SOS.get_SOS_kurtosis_estimates(kur_mat, deg_range) 
    f₄ₘᵢₙ_gs_df = gs.batch_min(f₄, n_stocks, grid_range)
    return Dict("SOS_estimates"=> f₄ₘᵢₙ_pop_df,
                "grid_search_estimates"=> f₄ₘᵢₙ_gs_df,
                "convex_quadratic_estimate"=> f₄ₘᵢₙ_cqr)
end

function get_𝒻_tilde(λ, mom_mats, f₁ₘₐₓ, f₂ₘᵢₙ, f₃ₘₐₓ, f₄ₘᵢₙ, grid_range)
    n_stocks = length(mom_mats["means_vector"])
    fₒₚₜ = get_fₒₚₜ(f₁ₘₐₓ, f₂ₘᵢₙ, f₃ₘₐₓ, f₄ₘᵢₙ)
    𝒻_tilde = obf.get_𝒻(λ, fₒₚₜ, mom_mats)
    𝒻_tilde_df = grid_search.batch_min(𝒻_tilde, n_stocks, grid_range)
    return Dict("grid_search_estimates"=> 𝒻_tilde_df)
end 

function get_fₒₚₜ(f₁ₘₐₓ, f₂ₘᵢₙ, f₃ₘₐₓ, f₄ₘᵢₙ)
    return (f₁ₘₐₓ["f₁ₘₐₓ_val"],
            f₂ₘᵢₙ["f₂ₘᵢₙ_val"],
            f₃ₘₐₓ["SOS_estimates"].obj_val[end],
            f₄ₘᵢₙ["SOS_estimates"].obj_val[end])
end

###################################################################################
function save_MVSK(MVSK_data, MVSK_bounds; delim=',')
    save_dir = get_save_dir(MVSK_data.n_stocks)
    mkdir(save_dir)
    save_MVSK_data(MVSK_data, save_dir)
    save_MVSK_bounds(MVSK_bounds, save_dir, delim=delim)
    return save_dir
end

function save_MVSK_data(MVSK_data, save_dir)
    data_dir = save_dir*"Data\\"
    mkdir(data_dir )
    for mom_mat ∈ ["means_vector", "kurtosis_matrix", "covariance_matrix", "skewness_matrix"]
        CSV.write(data_dir*mom_mat*".csv", DataFrame(MVSK_data.moment_matrices[mom_mat]),header = false)
    end
end

function save_MVSK_parameters(MVSK_bounds, save_dir; delim=',')
    error("Not coded yet")
end

function save_MVSK_bounds(MVSK_bounds, save_dir; delim=',')
    bounds_dir = save_dir*"Bounds\\"
    mkdir(bounds_dir)

    f₃ₘₐₓ_gs = MVSK_bounds.f₃ₘₐₓ["grid_search_estimates"]
    f₄ₘᵢₙ_gs = MVSK_bounds.f₄ₘᵢₙ["grid_search_estimates"]

    f₃ₘₐₓ_pop = MVSK_bounds.f₃ₘₐₓ["SOS_estimates"]
    f₄ₘᵢₙ_pop = MVSK_bounds.f₄ₘᵢₙ["SOS_estimates"]

    𝒻_tilde_gs = MVSK_bounds.𝒻_tilde["grid_search_estimates"]

    CSV.write(bounds_dir*"\\f₃ₘₐₓ_gs.csv", f₃ₘₐₓ_gs, header = true, delim=delim)
    CSV.write(bounds_dir*"\\f₃ₘₐₓ_pop.csv", f₃ₘₐₓ_pop, header = true, delim=delim)

    CSV.write(bounds_dir*"\\f₄ₘᵢₙ_gs.csv", f₄ₘᵢₙ_gs, header = true, delim=delim)
    CSV.write(bounds_dir*"\\f₄ₘᵢₙ_pop.csv", f₄ₘᵢₙ_pop, header = true, delim=delim)

    CSV.write(bounds_dir*"\\𝒻_tilde_gs.csv", 𝒻_tilde_gs, header = true, delim=delim)
end

get_time_now() = replace(string(now()), ":" => "-")
get_save_dir(n_stocks) = "assets\\"*get_time_now()*"n_stocks_$(n_stocks)\\"


end
