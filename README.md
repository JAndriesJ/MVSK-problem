# Mean-Variance-Skewness-Kurtosis (MVSK) optimization 
![[MVSK.png]]

## Introduction (Math and theory)


---
## How to install
- Download
- Install
- Instantiate
- Test installation

---

## Getting started
```Julia
include("MVSK.jl")            ; using .MVSK     

#-------------------------------------------- Data -------------------------------------------- 
pd = MVSK.get_processed_stock_data()
#--------------------------------------------single calc-------------------------------------------- 
λ = [0.2878  0.302  0.2831  0.1271]

simp_F_λ_opt = MVSK.get_F_λ_optimal_(pd.M,pd.V,pd.S,pd.K,λ)
simp_F_λ_opt.domain
simp_F_λ_opt.optimizer 
simp_F_λ_opt.opt_val
simp_F_λ_opt.opt_stat
simp_F_λ_opt.sol_time
simp_F_λ_opt.λ

box_F_λ_opt = MVSK.get_F_λ_optimal_(pd.M,pd.V,pd.S,pd.K,λ,1)
box_F_λ_opt.domain
box_F_λ_opt.optimizer 
box_F_λ_opt.opt_val
box_F_λ_opt.opt_stat
box_F_λ_opt.sol_time
box_F_λ_opt.λ

# Sparse one ??


#-------------------------------------------- -------------------------------------------- 
# MVSK.get_F_Λ_Pareto(save_path="default.csv",mesh_fineness=10,box_size=0,sub=0)
load_path = "C:\\Users\\jandr\\code_projects\\MVSK\\assets\\"* ["pareto_071122_40_simp_sparse_5.csv",
                                                                "pareto_261022_40_simp.csv",
                                                                "pareto_271022_40_box.csv"][2]
MVSK.plot_Pareto(load_path, sel=[4])
```



### How to use your own data.
All you need to do is replace the "...\MVSK\assets\stock_prices.csv" with another CSV that is ove the form

|date|Stock1|Stock2|Stock3|Stock4|...|
|---|---|---|---|---|---|
|1989-12-29|0.117203|0.352438|3.9375|3.48607 |1.752478|2.365775|1.766756|
|1990-01-02|0.123853|0.364733|4.125 |3.660858|1.766686|2.398184|1.766756|
|1990-01-03|0.124684|0.36405 |4.0   |3.660858|1.780897|2.356   |0.173216|
$\vdots$



## Technologies used (libraries & versions)


### Acknowledgments
This code was written as part of my [POEMA](http://poema-network.eu/) secondment at [NAG](https://www.nag.com/)

This project was funded in part by the Europeans Union's EU Framework Programme for Research and Innovation Horizon 2020 under Marie Skodowska-Curie Actions grant agreement 813211 (POEMA).

I would like to thank:
    1. my promoter Monique Laurent for her guidance and patience throughout the duration of my PhD.
    2. my contacts at NAG: Jan and Shuan for their support and guidance.










