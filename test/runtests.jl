using MVSK
using Test, JuMP, DataFrames

include("src\\stock_data.jl")
include("src\\moments.jl")
include("src\\constraints.jl")
include("src\\SDPmodel.jl")
include("src\\SDPoptimized.jl")
using .stock_data
using .moments 
using .constraints 
using .SDPmodel 
using .SDPoptimized

## stock_data

@testset "read_stock_data" begin
    read_stock_data = stock_data.read_data()
    @test typeof(read_stock_data) == DataFrames.DataFrame
    @test size(read_stock_data) == (501,20)
    @test names(read_stock_data)[2:6] == ["AAPL", "FB", "BABA", "AMZN", "GE"]
end

@testset "load_relative_returns_data_matrix" begin
    relative_returns_data_matrix_Google_Apple = stock_data.load_relative_returns_data_matrix(["GOOG","AAPL"])
    @test size(relative_returns_data_matrix_Google_Apple) == (500,2) 
    @test relative_returns_data_matrix_Google_Apple[37,:] == [0.016228290812141218, -0.0009086857825780688]


    relative_returns_data_matrix_first_six = stock_data.load_relative_returns_data_matrix([1:6 ...])
    @test size(relative_returns_data_matrix_first_six) == (500,6) 
    @test relative_returns_data_matrix_first_six[37,:] == [0.016228290812141218, -0.0009086857825780688, 0.0053498385640313, -0.002569334664405638, 0.0040069984249452816, 0.005640331169916136]
end

@testset "get_data_moment" begin
    relative_returns_data_matrix_first_six = stock_data.load_relative_returns_data_matrix([1:6 ...])
    moment_expo_2 = [0,0,0,0,1,1]
    moment_expo_3 = [1,0,0,0,1,1]
    moment_expo_4 = [3,0,0,0,0,1]
    @test stock_data.get_data_moment(relative_returns_data_matrix_first_six, moment_expo_2) == 3.74186338484832e-5
    @test stock_data.get_data_moment(relative_returns_data_matrix_first_six, moment_expo_3) == -4.150187958278357e-6
    @test stock_data.get_data_moment(relative_returns_data_matrix_first_six, moment_expo_4) == 1.477090175511437e-7
end

## moments

@testset "eᵢ" begin
    n = 7
    for k in 1:n
        e = moments.eᵢ(n,k)
        @test length(e) == n
        @test maximum(e) == 1
        @test minimum(e) == 0
        @test sum(e) == 1
    end
end

@testset "make_mom_expo_mat" begin
    mon_expo_mat = moments.make_mon_expo(3,(1,1);isle = true)
    @test [mon_expo_mat[i,i] for i in 1:4]  == [ [0, 0, 0],[2, 0, 0],[0, 2, 0],[0, 0, 2]]
    n,t = 4,2
    M_vec =  moments.make_mon_expo(n,t)
    @test length(M_vec) == binomial(n+t,t)
    M_mat =  moments.make_mon_expo(n,(t,t))
    @test M_vec == M_mat[1,:]
end

@testset "(Lx,α) -> var[x^α]" begin
    N,t = 5,4
    model = Model()
    moment_vec = moments.make_mon_expo(N,2*t)
    @variable(model, Lx[moment_vec]) # Create variables
    α_list  = [[0, 1, 0, 0, 0], [0, 0, 0, 4, 4], [0, 0, 0, 0, 8]]
    @test moments.get_Lxᵅ(Lx,α_list) == [Lx[[0, 1, 0, 0, 0]], Lx[[0, 0, 0, 4, 4]], Lx[[0, 0, 0, 0, 8]]]
    β_list  =  [[0, 6, 0, 1, 0],[3, 3, 0, 0, 2],[0, 2, 6, 0, 0],[2, 0, 0, 1, 3]]
    @test moments.get_Lxᵅ(Lx,α_list) == [Lx[[0, 6, 0, 1, 0]], Lx[[3, 3, 0, 0, 2]], Lx[[0, 2, 6, 0, 0]], Lx[[2, 0, 0, 1, 3]]]
end

## constraints
@testset "all the constraints" begin
    model = JuMP.Model()
    N,t = 5,3
    @variable(model, Lx[moments.make_mon_expo(N,2*t)]) # Create variables

    @test Lx[moments.make_mon_expo(N,0)[1,1]] == Lx[[0, 0, 0, 0, 0]]   ## L(1) 
    PSD_moment_matrix     = constraints.make_PSD_constraint(N,t,Lx) ## L([x]≦ₜ[x]ᵀ≦ₜ)
    @test size(PSD_moment_matrix) == (binomial(N+t,t),binomial(N+t,t))
    @test PSD_moment_matrix[end,end] == Lx[moments.eᵢ(N,N)*2*t]

    localzing_constraints = [con for con in constraints.make_localzing_constraints(N,t,Lx)] ## L(xᵢ[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0 ∀ i ∈ [N]
    @test length(localzing_constraints) == N
    @test size(localzing_constraints[1]) == (binomial(N+t-1,t-1),binomial(N+t-1,t-1))
    @test localzing_constraints[2][end,end] == Lx[[0, 1, 0, 0, 4]]

    loc_ideal_con         = constraints.make_localizing_ideal_constraint(N,t,Lx) ## L((1 - ∑ᴺᵢ₌₁xᵢ)[x]≦ₜ₋₁[x]ᵀ≦ₜ₋₁) ⪰ 0
    @test size(loc_ideal_con) == (binomial(N+t-1,t-1),binomial(N+t-1,t-1))
    @test loc_ideal_con[end,end] == Lx[[0, 0, 0, 0, 4]] - Lx[[1, 0, 0, 0, 4]] - Lx[[0, 1, 0, 0, 4]] - Lx[[0, 0, 1, 0, 4]] - Lx[[0, 0, 0, 1, 4]] - Lx[[0, 0, 0, 0, 5]]

    φ⁽³⁾ = constraints.make_objective_function(N,3,Lx) ## wᵀΦ⁽³⁾(w⊗w)
    φ⁽⁴⁾ = constraints.make_objective_function(N,4,Lx) ## wᵀΦ⁽⁴⁾(w⊗w⊗w)
    @test φ⁽³⁾.terms.keys == moments.get_Lxᵅ(Lx,moments.make_mon_expo(N,3,isle=false) )
    @test φ⁽⁴⁾.terms.keys == moments.get_Lxᵅ(Lx,moments.make_mon_expo(N,4,isle=false) )
end

## model
@testset "SDPmodel" begin
    N,t,k = 5,3,3
    SDP_model = SDPmodel.get_SDP_model(N,t,k)
end

@testset "SDPoptimized" begin
    SDP_model_opt = SDPoptimized.optimize_SDP(SDP_model)
    @test string(primal_status(SDP_model_opt)) == "FEASIBLE_POINT"
    @test string(dual_status(SDP_model_opt)) == "FEASIBLE_POINT"
    @test objective_value(SDP_model_opt) == 7.076637829729638e-19
end


@testset "MVSK.jl" begin
    # Write your tests here.




end
