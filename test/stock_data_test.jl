module stock_data_test
using Test
using DataFrames
using LinearAlgebra

include(pwd()*"/src/stock_data.jl")            ; using .stock_data   
# fieldnames(typeof(pd)) == (:R, :M, :std, :V, :S, :S_mat, :K, :K_mat, :R_max, :R_min, :R_max_max, :R_min_min, :n, :m)
# df = stock_data.read_stock_csv()
function run_tests()
    @testset "read_stock_csv" begin
        df_raw = stock_data.read_stock_csv()
        @test typeof(df_raw) == DataFrame
        @test size(df_raw) == (501,20)
    end

    @testset "calc_relative_returns" begin
        R = stock_data.calc_relative_returns(df_raw)
        @test size(R) == (500,20)
        @test R[1] == 0.010026330698287287
    end

    @testset "mean and std" begin
        M   = get_means_vector(R) 
        length(M) == 20
        std = get_standard_deviation(R) 
        @test all(std .≥ 0)
    end

    @testset "sizeand ent_std_prices" begin
        m,n = size(R)
        @test m == 500
        @test n == 20

        R     = calc_cent_std_prices(R,M,std,m,n) 
        @test size(R) ==  (500,20)
    end

    @testset "variance" begin
        V     = calc_V(R,m)
        @test isapprox(V, (1/m)*sum([R[i,:]*R[i,:]' for i ∈ 1:m]) ,atol=1e-8)
        @test issymmetric(V) 
        @test minimum(eigvals(V)) ≥ -1e-12
    end

    @testset "skewness" begin    
        S_mat = calc_S_mat(R,m)
        S     = calc_S(R,m,n)
        @test S_mat[1] == S[1]
        @test size(S_mat) == (n^2, n)
        @test size(S) == (n, n, n)
    end

    @testset "kurtosis" begin
        K_mat = calc_K_mat(R,m) 
        K     = calc_K(R,m,n)
        @test K_mat[1] == K[1]
        @test issymmetric(K_mat)
        @test minimum(eigvals(K_mat)) ≥ -1e-12
    end

    @testset "maxima" begin    
        R_max = maximum(R,dims=2) # [maximum(R[t,:]) for t ∈ 1:m]   
        R_min = minimum(R,dims=2) 
        R_max_max = maximum(R_max)
        R_min_min = minimum(R_min)

        @test all(R_min .≤ R_max)
        @test maximum(R_max) == R_max_max 
        @test R_min_min == minimum(R_min)
    end 

    # read_proc_data(proc_data_path = data_dir*"proc_data") =  deserialize(proc_data_path)
    # write_sproc_data(pd::proc_data; proc_data_path = data_dir*"proc_data") =  serialize( proc_data_path, pd)

end



end










