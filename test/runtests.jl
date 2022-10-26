module test

include("stock_data_test.jl");  using .stock_data_test   
include("pgd_test.jl")       ;  using .pgd_test
include("spaces_test.jl")    ;  using .spaces_test
include("λ_char_test.jl")    ;  using .λ_char_test

include("obj_test.jl")       ; using .obj_test
include("batch_test.jl")     ; using .batch_test



stock_data_test.run_tests()
pgd_test.run_tests()
spaces_test.run_tests()
λ_char_test.run_tests()
obj_test.run_tests()
batch_test.run_tests()

end
# using MVSK
# using Test, JuMP, DataFrames

# include("src\\stock_data.jl")
# include("src\\moments.jl")
# include("src\\constraints.jl")
# include("src\\SDPmodel.jl")
# include("src\\SDPoptimized.jl")
# using .stock_data
# using .moments 
# using .constraints 
# using .SDPmodel 
# using .SDPoptimized

# ## stock_data

# @testset "read_stock_data" begin
#     read_stock_data = stock_data.read_data()
#     @test typeof(read_stock_data) == DataFrames.DataFrame
#     @test size(read_stock_data) == (501,20)
#     @test names(read_stock_data)[2:6] == ["AAPL", "FB", "BABA", "AMZN", "GE"]
# end

# @testset "load_relative_returns_data_matrix" begin
#     relative_returns_data_matrix_Google_Apple = stock_data.load_relative_returns_data_matrix(["GOOG","AAPL"])
#     @test size(relative_returns_data_matrix_Google_Apple) == (500,2) 
#     @test relative_returns_data_matrix_Google_Apple[37,:] == [0.016228290812141218, -0.0009086857825780688]


#     relative_returns_data_matrix_first_six = stock_data.load_relative_returns_data_matrix([1:6 ...])
#     @test size(relative_returns_data_matrix_first_six) == (500,6) 
#     @test relative_returns_data_matrix_first_six[37,:] == [0.016228290812141218, -0.0009086857825780688, 0.0053498385640313, -0.002569334664405638, 0.0040069984249452816, 0.005640331169916136]
# end

# @testset "get_data_moment" begin
#     relative_returns_data_matrix_first_six = stock_data.load_relative_returns_data_matrix([1:6 ...])
#     moment_expo_2 = [0,0,0,0,1,1]
#     moment_expo_3 = [1,0,0,0,1,1]
#     moment_expo_4 = [3,0,0,0,0,1]
#     @test stock_data.get_data_moment(relative_returns_data_matrix_first_six, moment_expo_2) == 3.74186338484832e-5
#     @test stock_data.get_data_moment(relative_returns_data_matrix_first_six, moment_expo_3) == -4.150187958278357e-6
#     @test stock_data.get_data_moment(relative_returns_data_matrix_first_six, moment_expo_4) == 1.477090175511437e-7
# end





# ## model
# @testset "SDPmodel" begin
#     N,t,k = 5,3,3
#     SDP_model = SDPmodel.get_SDP_model(N,t,k)
# end

# @testset "SDPoptimized" begin
#     SDP_model_opt = SDPoptimized.optimize_SDP(SDP_model)
#     @test string(primal_status(SDP_model_opt)) == "FEASIBLE_POINT"
#     @test string(dual_status(SDP_model_opt)) == "FEASIBLE_POINT"
#     @test objective_value(SDP_model_opt) == 7.076637829729638e-19
# end


# @testset "MVSK.jl" begin
#     # Write your tests here.
# end
