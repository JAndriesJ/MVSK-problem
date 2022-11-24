module test

include("stock_data_test.jl");  using .stock_data_test   
include("spaces_test.jl")    ;  using .spaces_test
include("λ_char_test.jl")    ;  using .λ_char_test
include("obj_test.jl")       ; using .obj_test
# include("pareto_set_test.jl"); using .pareto_set_test

stock_data_test.run_tests()
spaces_test.run_tests()
λ_char_test.run_tests() 
obj_test.run_tests()
# pareto_set_test.run_tests()

end
