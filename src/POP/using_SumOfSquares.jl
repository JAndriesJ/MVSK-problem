module using_SumOfSquares

using DynamicPolynomials
using SumOfSquares
using JuMP
using LinearAlgebra
using Test
import CSDP

export get_SOS_bound,
       run_tests  

"""
max  γ 
s.t. xᵀ Ψ (x ⊗ x) or (x ⊗ x)ᵀ Ψ (x ⊗ x) - γ ∈ SOS
x ∈ Δⁿ
"""
function get_SOS_bound(matrix,k)  #
    n,m = size(matrix)
    @assert m == n^2 || m == n
    n == m ? m = n = Int(sqrt(n)) : 1
    @assert typeof(n) <: Int ; @assert 1 < n < 11 ; 
    @polyvar x[1:n] 
    f =  m == n^2 ? x[1:n]'*matrix*kron(x[1:n],x[1:n]) : kron(x[1:n],x[1:n])'*matrix*kron(x[1:n],x[1:n]) 
    ## This is a bootleg solution
    if n== 2
        S = @set x[1] >= 0 && x[2] >= 0 && (x[1:n]'*ones(n,1))[1] >= 1
    elseif n== 3
        S = @set x[1] >= 0 && x[2] >= 0 && x[3] >= 0 && (x[1:n]'*ones(n,1))[1] >= 1
    elseif n== 4
        S = @set x[1] >= 0 && x[2] >= 0 && x[3] >= 0 && x[4] >= 0 && (x[1:n]'*ones(n,1))[1] >= 1
    elseif n== 5
        S = @set x[1] >= 0 && x[2] >= 0 && x[3] >= 0 && x[4] >= 0 && x[5] >= 0 && (x[1:n]'*ones(n,1))[1] >= 1
    elseif n == 6
        S = @set x[1] >= 0 && x[2] >= 0 && x[3] >= 0 && x[4] >= 0 && x[5] >= 0 && x[6] >= 0 &&  (x[1:n]'*ones(n,1))[1] >= 1    
    elseif n == 7
        S = @set x[1] >= 0 && x[2] >= 0 && x[3] >= 0 && x[4] >= 0 && x[5] >= 0 && x[6] >= 0 && x[7] >= 0 &&  (x[1:n]'*ones(n,1))[1] >= 1
    elseif n == 8
        S = @set x[1] >= 0 && x[2] >= 0 && x[3] >= 0 && x[4] >= 0 && x[5] >= 0 && x[6] >= 0 && x[7] >= 0 && x[8] >= 0 && (x[1:n]'*ones(n,1))[1] >= 1
    elseif n == 9
        S = @set x[1] >= 0 && x[2] >= 0 && x[3] >= 0 && x[4] >= 0 && x[5] >= 0 && x[6] >= 0 && x[7] >= 0 && x[8] >= 0 && x[9] >= 0 && (x[1:n]'*ones(n,1))[1] >= 1    
    elseif n == 10
        S = @set x[1] >= 0 && x[2] >= 0 && x[3] >= 0 && x[4] >= 0 && x[5] >= 0 && x[6] >= 0 && x[7] >= 0 && x[8] >= 0 && x[9] >= 0 && x[10] >= 0 && (x[1:n]'*ones(n,1))[1] >= 1
    end

    solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
    model = SOSModel(solver)
    @variable(model, α)
    @objective(model, Max, α)
    @constraint(model, c, f >= α, domain = S, maxdegree = k)
    optimize!(model)
    @show termination_status(model)
    return objective_value(model), primal_status(model), dual_status(model), solve_time(model)

end

function run_tests()
    @testset "get_SOS_bound skewness for random data" begin
        for n ∈ 2:10
            skewness_matrix = rand(n,n^2) 
            get_SOS_bound(skewness_matrix,3)
        end
    end

    @testset "get_SOS_bound kurtosis for random data" begin
        for n ∈ 2:10
            n = 4
            R =  rand(n,n^2)
            kurtosis_matrix = R'*R
            sos_bound = get_SOS_bound(kurtosis_matrix,4)
            @test (1/16)*sos_bound ≤ minimum(diag(kurtosis_matrix)) # 
        end
    end

    @testset "special cases" begin
        skewness_matrix = kron(ones(1,4),I(4))
        @test get_SOS_bound(skewness_matrix,3) ≈ 0.25

        kurtosis_matrix = kron(I(4),I(4))
        @test get_SOS_bound(kurtosis_matrix,4) ≈ 0.06249999986000854
    end
end

end
"""
Sanity check:
Theory: 
1. Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.
2. Lasserre, J. B. Moments, positive polynomials and their applications. World Scientific, 2009.
3. Laurent, M. Sums of squares, moment matrices and optimization over polynomials Emerging applications of algebraic geometry, Springer, 2009, 157-270.

## Example https://jump.dev/SumOfSquares.jl/stable/generated/Polynomial%20Optimization/polynomial_optimization/

min x³ - x² + 2xy - y² + y³ 
s.t. x ≥ 0, y ≥ 0, x + y ≥ 1
Global optimals: (1,0), (0,1)
Local optimal: (1/2,1/2)

max γ
s.t. x³ - x² + 2xy - y² + y³ - γ ∈ SOS(x ≥ 0, y ≥ 0, x + y ≥ 1) 

```
@polyvar x y
p = x^3 - x^2 + 2x*y -y^2 + y^3
S = @set x >= 0 && y >= 0 && x + y >= 1
p(1,0), p(1//2,1//2), p(0,1)

solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver)
@variable(model, γ)
@objective(model, Max, γ)
@constraint(model, c5, p >= γ, domain = S, maxdegree = 5)
optimize!(model)
@show termination_status(model)
@show objective_value(model)

ν5 = moment_matrix(c5)
ν5.Q
ν5.basis.monomials
extractatoms(ν5, 1e-3)
```
"""




