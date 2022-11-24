module stock_data_test
using Test
using DataFrames
using LinearAlgebra
 
include(pwd()*"/src/stock_data.jl")            ; using .stock_data   

function run_tests()
    @testset "read_stock_csv" begin
        global df_raw = stock_data.read_stock_csv()
        @test typeof(df_raw) == DataFrame
        @test size(df_raw) == (501,20)
    end

    @testset "calc_relative_returns" begin
        global R = stock_data.calc_relative_returns(df_raw)
        @test size(R) == (500,20)
        @test R[1] == 0.010026330698287287
    end

    @testset "mean and std" begin
        global M = stock_data.get_means_vector(R) 
        length(M) == 20
        global std = stock_data.get_standard_deviation(R) 
        @test all(std .≥ 0)
    end

    @testset "cent_std_prices" begin
        global m = 500
        global n = 20
        Rcs = stock_data.calc_cent_std_prices(R,M,std) 
        @test isapprox(stock_data.get_means_vector(Rcs), zeros(n), atol=1e-14)
        @test isapprox(stock_data.get_standard_deviation(Rcs), ones(n), atol=1e-14)
        @test size(Rcs) ==  (500,20)
    end

    @testset "variance" begin
        V = stock_data.calc_V(R,m)
        isapprox(std, sqrt.(diag(V)) ,atol=1e-3)
        @test isapprox(V, (1/m)*sum([R[i,:]*R[i,:]' for i ∈ 1:m]) ,atol=1e-14)
        @test LinearAlgebra.issymmetric(V) 
        @test minimum(eigvals(V)) ≥ -1e-12
    end

    @testset "skewness" begin    
        S_mat = stock_data.calc_S_mat(R,m)
        S     = stock_data.calc_S(R,m,n)
        @test S_mat[1] ≈ S[1]
        @test size(S_mat) == (n^2, n)
        @test size(S) == (n, n, n)
    end

    @testset "kurtosis" begin
        K_mat = stock_data.calc_K_mat(R,m) 
        K     = stock_data.calc_K(R,m,n)
        @test K_mat[1] ≈ K[1]
        @test LinearAlgebra.issymmetric(K_mat)
        @test minimum(eigvals(K_mat)) ≥ -1e-12
    end

    @testset "maxima" begin    
        R_max = maximum(R,dims=2)   
        R_min = minimum(R,dims=2) 
        R_max_max = maximum(R_max)
        R_min_min = minimum(R_min)

        @test all(R_min .≤ R_max)
        @test maximum(R_max) == R_max_max 
        @test R_min_min == minimum(R_min)
    end 

    @testset "proc_data" begin
        pd = stock_data.proc_data(R)
        @test fieldnames(typeof(pd)) == (:R, :M, :std, :V, :S, :S_mat, :K, :K_mat, :R_max, :R_min, :R_max_max, :R_min_min, :n, :m)
    end

end



end










