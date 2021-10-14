# Mean-Variance-Skewness-Kurtosis (MVSK) optimization 
![[MVSK.png]]

## Introduction 
This repository provides polynomial optimization tools for the solving of the MVSK problem from portfolio optimization.
It follows the general paradigm of moment based hierarchies of approximations for polynomial optimization problems.
This code was written as part of my [POEMA](http://poema-network.eu/) secondment at [NAG](https://www.nag.com/)
For details of the mathematics see [arXiv paper](url)

## How to install
The code should be stand alone:
git clone ....
Pkg.instantiate()
ready to roll


## Background:
The sole purpose of this code is the computation of the the following polynomial optimization problem (POP)
$$\begin{equation} \tag{3}
\begin{split} 
\psi^{(3)}_{\pm}= \max \textbackslash \min & ~~ w^T\Phi^{(3)}(w\otimes w) \\
s.t. &~~  w \in \Delta^N \\
\end{split}
 \end{equation}$$

$$\begin{equation} \tag{4}
\begin{split} 
\psi^{(4)}_{\pm}= \max \textbackslash \min & ~~ w^T\Phi^{(4)}(w\otimes w \otimes w) \\
s.t. &~~  w \in \Delta^N \\
\end{split}
 \end{equation}$$
Since POPs are NP hard to solve in general, even on the simplex, we resort to approximations. We use the moment based hierarchy of approximations:
$$
\begin{equation} \tag{mom3}
\begin{split} 
\min & ~~  L(f) \\
s.t. &~~ L(1) = 1\\
&~~ M(g_0L) := L([\bw][\bw]^T) \succeq 0 \\
&~~ M(g_iL) := L(x_i [\bw][\bw]^T) \succeq 0 ~\forall~i \in [N]\\
&~~ M(hL):= L( (1 - \sum_{i\in [N]} x_i) [\bw][\bw]^T) = 0\\
\end{split}
\end{equation}
$$
Where the following are user set parameters:
- $f(\bw) =  w^T\Phi^{(k)}(w\otimes ... \otimes w)$
    - $\Phi^{(k)}$ is the $k^{th}$ moment matrix.
- $N$ is the number of stocks.
- $t$ is the level of the hierarchy.

## Getting started
```Julia
## Running a single model:
Number_of_stocks = 10
Level_of_Hierarchy = 2 # must be higher than 1
k = 3 # 3 for skewnes, 4 for kurtosis
mod     = SDPmodel.get_SDP_model(Number_of_stocks,Level_of_Hierarchy,3)
mod_opt = SDPoptimized.optimize_SDP(mod)


Primal_status    = string(primal_status(mod_opt)) 
Dual_status      = string(dual_status(mod_opt))
objective_value  = JuMP.objective_value(mod_opt)
computation_time = JuMP.solve_time(mod_opt)
```

```Julia
## Running the model for a batch of parameters:
## Skewness
Number_of_stocks_list = [5:15 ...]
Level_of_Hierarchy_list = [2:3 ...] # must be higher than 1
k = 3 # 3 for skewnes, 4 for kurtosis
SDPoptimized.batch_optimize_SDP(Number_of_stocks_list,Level_of_Hierarchy_list,k)
## Kurtosis
Number_of_stocks_list = [5:15 ...]
Level_of_Hierarchy_list = [2:3 ...] # must be higher than 1
k = 4 # 3 for skewnes, 4 for kurtosis
SDPoptimized.batch_optimize_SDP(Number_of_stocks_list,Level_of_Hierarchy_list,k)
```



### How to get your own data.
All you need to do is replace the "...\MVSK\assets\stock_prices.csv" with another CSV that is ove the form

|date|Stock1|Stock2|Stock3|Stock4|...|
|---|---|---|---|---|---|
|1989-12-29|0.117203|0.352438|3.9375|3.48607 |1.752478|2.365775|1.766756|
|1990-01-02|0.123853|0.364733|4.125 |3.660858|1.766686|2.398184|1.766756|
|1990-01-03|0.124684|0.36405 |4.0   |3.660858|1.780897|2.356   |0.173216|
$\vdots$



## Technologies used (libraries & versions, helps recruiters)


### Acknowledgments
This project was funded in part by the Europeans Union's EU Framework Programme for Research and Innovation Horizon 2020 under Marie Skodowska-Curie Actions grant agreement 813211 (POEMA).

I would like to thank:
    1. my promoter Monique Laurent for her guidance and patience throughout the duration of my PhD.
    2. my contacts at NAG: Jan and Shuan for their support and guidance.










