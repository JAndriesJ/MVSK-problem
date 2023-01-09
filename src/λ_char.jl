module λ_char

export  λ₃_bound,
        λ₃_bound_Δ

"""max_λ₃ s.t. 2λ₂ - 6λ₃*y +12λ₄*y^2 ≥ 0 for all y"""
λ₃_bound(λ₂, λ₄)     = sqrt((8/3)*λ₂*λ₄)
λ₃_bound(λ₂, λ₃, λ₄) = λ₃_bound(λ₂, λ₄) ≥ λ₃
λ₃_bound(λ)          = λ₃_bound(λ[2], λ[4]) ≥ λ[3]

"""max_λ₃ s.t. 2λ₂ - 6λ₃*y +12λ₄*y^2 ≥ 0 for all y ...simplex or box [-1,1]^n"""
λ₃_bound_Δ(λ₂, λ₄, R_up::Real)  =  (sqrt(λ₂/(6*λ₄)) > R_up) ? (λ₂ + 6*λ₄*(R_up^2))/(3*R_up) : λ₃_bound(λ₂, λ₄)
λ₃_bound_Δ(λ₂, λ₃, λ₄, R_up::Real) =  λ₃_bound_Δ(λ₂, λ₄, R_up) ≥ λ₃
λ₃_bound_Δ(λ, R_up::Real) =  λ₃_bound_Δ(λ[2], λ[4], R_up) ≥ λ[3]




end



