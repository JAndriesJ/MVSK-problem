module stock_data

using Statistics
using CSV, DataFrames
using LinearAlgebra ; const la = LinearAlgebra 
using Serialization

export  read_stock_csv,
        read_proc_data,
        write_proc_data,
        get_proc_data

global data_dir = pwd()*"\\assets\\"

function get_proc_data(csv_path = data_dir*"stock_prices.csv")
    df = read_stock_csv(csv_path)
    R  = calc_relative_returns(df) 
    pd = process_data(R)
    # write_proc_data(pd) 
    return pd
end

struct proc_data
    R 
    M 
    std
    V
    S
    S_mat
    K
    K_mat
    R_max
    R_min
    R_max_max
    R_min_min
    n
    m
end

function process_data(R::Matrix{Float64}) 
    M   = get_means_vector(R) 
    std = get_standard_deviation(R) 
    m,n = size(R)

    R     = calc_cent_std_prices(R,M,std)  
    V     = calc_V(R,m)
    S_mat = calc_S_mat(R,m)
    S     = calc_S(R,m,n)
    K_mat = calc_K_mat(R,m) 
    K     = calc_K(R,m,n)

    R_max = maximum(R,dims=2) # [maximum(R[t,:]) for t ∈ 1:m]   
    R_min = minimum(R,dims=2) 
    R_max_max = maximum(R_max)
    R_min_min = minimum(R_min)

    pd =  proc_data(R,
                    M,
                    std,
                    V,
                    S,
                    S_mat,
                    K,
                    K_mat,
                    R_max,
                    R_min,
                    R_max_max,
                    R_min_min,
                    n,
                    m)

    return pd
end
get_means_vector(R) = vec(Statistics.mean(R, dims=1))
get_covariance_matrix(R) = Statistics.cov(R, dims=1)
get_standard_deviation(R)  = sqrt.(diag(get_covariance_matrix(R)))
calc_cent_std_prices(R,M,std) =  (R - repeat(reshape(M,1,size(R)[2]),size(R)[1])) ./ repeat(std',size(R)[1])

calc_V(R,m) = (1/(m))*sum([R[i,:]*R[i,:]' for i ∈ 1:m])
calc_S_mat(R,m) = (1/m)*sum([la.kron(R[i,:],R[i,:])*R[i,:]' for i ∈ 1:m])
calc_S(R,m,n) = reshape([sum([R[x,i]*R[x,j]*R[x,k] for x ∈ 1:m])/m for i in 1:n for j in 1:n for k in 1:n],n,n,n)
calc_K_mat(R,m) = (1/m)*sum([la.kron(R[i,:],R[i,:])*la.kron(R[i,:],R[i,:])' for i ∈ 1:m])
calc_K(R,m,n) = reshape([sum([R[x,i]*R[x,j]*R[x,k]*R[x,l] for x ∈ 1:m])/m for i in 1:n for j in 1:n for k in 1:n for l in 1:n],n,n,n,n)

function calc_relative_returns(df::DataFrame) 
    R = Matrix(df)
    return diff(R,dims=1)./ R[1:end-1,:]
end

function read_stock_csv(csv_path = data_dir*"stock_prices.csv")
    df = CSV.read(csv_path, DataFrame)
    return select!(df[end-500:end,:], Not(:date))
end

read_proc_data(proc_data_path = data_dir*"proc_data") =  deserialize(proc_data_path)
write_proc_data(pd::proc_data; proc_data_path = data_dir*"proc_data") =  serialize( proc_data_path, pd)

end



