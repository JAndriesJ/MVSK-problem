module Data_moments
using LinearAlgebra, Test

export get_𝒻,
       get_f₁,
       get_f₂,
       get_f₃,
       get_f₄,
       get_means_vector,
       get_covariance_matrix,
       get_skewness_matrix,
       get_kurtosis_matrix,
       run_tests

function get_𝒻(λ,ϕ,R)
    f₁  = get_f₁(R); f₂  = get_f₂(R); f₃  = get_f₃(R); f₄  = get_f₄(R) 
    𝒻(x) =  f₃(x) ≤ ϕ[3] && f₄(x) ≥ ϕ[4] ?  x -> (1 - (f₁(x)/ϕ[1]))^λ[1] + ((f₂(x)/ϕ[2]) - 1)^λ[2] + ((f₃(x)/ϕ[3]) - 1)^λ[3] + (1 - (f₄(x)/ϕ[1]))^λ[4] : x -> 0
    return 𝒻
end       

get_f₁(R) = f(x) = (get_means_vector(R)*x)[1]
get_f₂(R) = f(x) = x'*get_covariance_matrix(R)*x
get_f₃(R) = f(x) = x'*get_skewness_matrix(R)*kron(x,x)
get_f₄(R) = f(x) = kron(x,x)'*get_kurtosis_matrix(R)*kron(x,x)

get_means_vector(R) = sum(R, dims=1)/(m(R)-1)
get_covariance_matrix(R) = [ sum(R[:,i].*R[:,j])  for i ∈ 1:n(R), j ∈ 1:n(R)] ./ (m(R)-1)^2
get_skewness_matrix(R)   = [ sum(R[:,i].*R[:,jk[1]].*R[:,jk[2]]) for i ∈ 1:n(R), jk ∈  [get_ij_pairs(R)...]  ]./ (m(R)-1)^3
get_kurtosis_matrix(R)   = [ sum(R[:,ij[1]].*R[:,ij[2]].*R[:,kl[1]].*R[:,kl[2]]) for ij ∈ [get_ij_pairs(R)...], kl ∈  [get_ij_pairs(R)...]  ] ./ (m(R)-1)^4

## utility
m(R) = size(R)[1] ; n(R) = size(R)[2] # R ∈ ℝᵐˣⁿ
get_ij_pairs(R) = [ (j,i)  for i ∈ 1:n(R), j ∈ 1:n(R)] #

# Tests
function run_tests()
    @testset "special case" begin
        M = 100
        N = 10
        R = ones(M,N)
        means_vector      = get_means_vector(R)
        (M/(M-1))*ones(N,1)
        @test means_vector == (M/(M-1))*ones(N,1)'
        covariance_matrix = get_covariance_matrix(R)
        @test covariance_matrix == (M/(M-1)^2)*ones(N,N)
        skewness_matrix = get_skewness_matrix(R)
        @test skewness_matrix == (M/(M-1)^3)*ones(N,N^2)
        kurtosis_matrix   = get_kurtosis_matrix(R)
        @test kurtosis_matrix == (M/(M-1)^4)*ones(N^2,N^2)
    end

    @testset "moment matrices sizes" begin
        M = 100
        for N ∈ rand(7:12,6)
            R = rand(M,N)
            means_vector      = get_means_vector(R)
            @test length(means_vector) == N
            covariance_matrix = get_covariance_matrix(R)
            @test size(covariance_matrix) == (N,N)
            skewness_matrix = get_skewness_matrix(R)
            @test size(skewness_matrix) == (N,N^2)
            kurtosis_matrix   = get_kurtosis_matrix(R)
            @test size(kurtosis_matrix) == (N^2,N^2)
        end
    end

    @testset "moment functions vs moment matrices" begin
        M = 100
        for N ∈ rand(7:12,6)
            R = rand(M,N)
            N
            M
            get_f₁(R)(ones(N))
            @test get_f₁(R)(ones(N)) ≈ sum(get_means_vector(R))
            @test get_f₂(R)(ones(N)) ≈ sum(get_covariance_matrix(R))
            @test get_f₃(R)(ones(N)) ≈ sum(get_skewness_matrix(R))
            @test get_f₄(R)(ones(N)) ≈ sum(get_kurtosis_matrix(R))
        end
    end
end

end


