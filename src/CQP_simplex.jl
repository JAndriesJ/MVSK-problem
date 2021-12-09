module CQP_simplex # Convex Quadratic programming on the simplex
using Ipopt, JuMP, LinearAlgebra, Test
export get_variance_min,
       run_tests  

"""
Solve the Convex Quadratic programming on the simplex:
    min xᵀQx
    s.t. x ∈ Δⁿ (the n-dimensional simplex)
    where Q ⪰ 0 is assumed but not requered.
"""     
function get_variance_min(covariance_matrix)
    n = size(covariance_matrix)[1]
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x[1:n] >= 0); @constraint(model, sum(x[1:n]) == 1) # x ∈ Δⁿ
    @objective(model, Min, x' * covariance_matrix * x )
    optimize!(model)
    return objective_value(model)
end

function get_local_min_𝒻(λ,f_para,n)
    return error
end


## Tests:
function run_tests()
    @testset "Ipopt on a linear problem" begin
        model = Model(Ipopt.Optimizer)
        @variable(model, x >= 1)
        @objective(model, Min, x + 0.5)
        optimize!(model)
        @test objective_value(model) ≈ 1.5
    end

    @testset "Ipopt on a degree 3 problem" begin    
        model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
        @variable(model, a >= 0); @variable(model, b >= 0)
        @constraint(model, a + b >= 1)
        @NLobjective(model, Min, a^3 - a^2 + 2a*b - b^2 + b^3)
        optimize!(model)
        @show termination_status(model)
        @show value(a)
        @show value(b)
        objective_value(model) ≈ 0.24999999500590336
    end    

    @testset "get_variance_min toy cases" begin
        @test get_variance_min([0,0,0]*[0,0,0]') == 0
        @test get_variance_min([1,1,1]*[1,1,1]') ≈ 1
        @test get_variance_min([1,0,0]*[1,0,0]') ≈ get_variance_min([0,1,0]*[0,1,0]') 
        @test get_variance_min([0,1,0]*[0,1,0]') ≈ get_variance_min([0,0,1]*[0,0,1]')
        @test get_variance_min([0,0,1]*[0,0,1]') ≈ 4.3652450324270914e-9
        @test get_variance_min([1,0,0]*[0,0,-1]') ≈ -0.25
    end

    @testset "get_variance_min_random variance" begin
        for n in rand(4:10,5)
            R_square = rand(n,n)
            R = R_square*R_square'
            min_diag = minimum([R[i,i] for i ∈ 1:n])
            @test get_variance_min(R) ≤ min_diag || get_variance_min(R) ≈ min_diag 
        end
    end
end

end
 

"""
Sanity check:
Theory: https://coin-or.github.io/Ipopt/
Ipopt (Interior Point Optimizer) is an open source software package for large-scale nonlinear optimization.
It can be used to solve general nonlinear programming problems of the form:
min f(x)
s.t. gᴸ ≤ g(x) ≤ gᵁ
     xᴸ ≤   x  ≤ xᵁ
The functions f(x) and g(x) can be linear or nonlinear and convex or non-convex (but should be twice continuously differentiable)
Note that equality constraints of the form gᵢ(x)=̄gᵢ can be specified by setting gᴸᵢ=gᵁᵢ=̄gᵢ.

Ipopt implements an ***interior point line search filter method that aims to find a local solution of (NLP)***.
The mathematical details of the algorithm can be found in several publications
1. J. Nocedal, A. Wächter, and R.A. Waltz. Adaptive barrier strategies for nonlinear interior methods. SIAM Journal on Optimization, 19(4):1674–1693, 2008. preprint at http://www.optimization-online.org/DB_HTML/2005/03/1089.html.
2. A. Wächter. An Interior Point Algorithm for Large-Scale Nonlinear Optimization with Applications in Process Engineering. PhD thesis, Carnegie Mellon University, Pittsburgh, PA, USA, January 2002.
3. A. Wächter and L.T. Biegler. On the implementation of a primal-dual interior point filter line search algorithm for large-scale nonlinear programming. Mathematical Programming, 106(1):25–57, 2006. preprint at http://www.optimization-online.org/DB_HTML/2004/03/836.html.
4. A. Wächter and L.T. Biegler. Line search filter methods for nonlinear programming: Motivation and global convergence. SIAM Journal on Optimization, 16(1):1–31, 2005.
5. A. Wächter and L.T. Biegler. Line search filter methods for nonlinear programming: Local convergence. SIAM Journal on Optimization, 16(1):32–48, 2005.

## Example see the tests above



"""