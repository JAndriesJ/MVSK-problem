using Pkg
Pkg.status()


using Plots ####### Pkg.add("Plots")
using Distributions # Pkg.add("Distributions")


N(x,μ=0,σ=1) = 1/(σ*√(2π)) * exp(-0.5*((x-μ)/σ)^2)


anim = Animation()
for σ in 0.5:0.1:4
    plot(-4:0.01:4,map(x->N(x,0,σ),-4:0.01:4),
    xticks=[],
     xlabel="μ",ylabel="Density",
     linewidth=3, label = "σ=$σ")
    frame(anim)
end
gif(anim,"Variance.gif",fps = 3)



Pkg.add("NumericalIntegration")
using NumericalIntegration
x = collect(-100 : 0.001 : α*(-))
y = sin.(x)

# integrate using the default Trapezoidal method


function sN(x,μ=0,σ=1,α=0)
    t = collect(-100 : 0.001 : α*((x-μ)/σ))
    t = collect(-100 : 0.001 : α*((x-μ)/σ))
    y =  1/(√(2π)) * exp.(-0.5.*t.^2)
    return N(x) * (integrate(t, y ) )
end



α = 0
plot(-4:0.01:4,map(x->sN(x,0,1,α),-4:0.01:4),
xticks=[],
 xlabel="μ",ylabel="Density",
 linewidth=3, label = "α=$α")


anim = Animation()
for α in -2:0.1:2
    plot(-4:0.01:4,map(x->sN(x,0,1,α),-4:0.01:4),
    xticks=[],
    xlabel="μ",ylabel="Density",
    linewidth=3, label = "α=$α")
    frame(anim)
end
gif(anim,"Skewness.gif",fps = 5)
 
