# Research into signomial optimization
## Date modified: 2021-11-21
### Objectives
1. Find the current best convergences result:
	1. That is the the theorem with the most general conditions for convergence of a hierarchy of SAGE-relaxations of a given signomial program.
	2. Known rates of convergence.
2. Create a list of applications of signomial optimization with associated references.
3. Find a minimal working example of numerical code for signomial optimization.

---
### Sources
#### Papers
- 1968,Duffin,Peterson,Zener, Geometric Programming Theory and Application.pdf ^DPZ68
- 1996,Sui,Wang,Second-order method of generalized geometric programming for spatial frame optimization ^SW96
- 1997,Maranas,Floudas,Global Optimization in Generalized Geometric Programming ^MF97
- 2005,Chiang,Geometric Programming For Communication Systems.pdf ^C05
- 2014,Chandrasekaran,Shah,Relative Entropy Relaxations for Signomial Optimization (ArXiv).pdf
- 2014,Chandrasekaran,Shah,Relative Entropy Relaxations for Signomial Optimization.pdf ^CS14
- 2014,Gzyl,Inverardi,Tagliani,Fractional Moments and Maximum Entropy Geometric Meaning.pdf ^GIT14
- 2014,Xu, Global optimization of signomial geometric programming problems.pdf ^X14
- 2016,Chandrasekaran,Shah,Relative entropy optimization and its applications.pdf ^16CS
- 2016,Iliman,de Wolf,Amoebas, nonnegative polynomials and sums of squares supported on circuits.pdf ^IW16
- 2017,Dressler,Iliman, De Wolff,an Approach to Constrained Polynomial Optimization.pdf ^DIW17
- 2019,Muller,Hofbauer,Regensburger,On the Bijectivity of Families of Exponential Generalized Polynomial Maps ^MHR19
- 2020,Katthan, Naumann, Theobald,a Unified Framework of Sage and Sonc Polynomials And.pdf ^KNT20
- 2020,MURRAY,NAUMANN,THEOBALD,SUBLINEAR CIRCUITS AND THE CONSTRAINED SIGNOMIAL NONNEGATIVITY PROBLEM.pdf ^MNT20
- 2020,Murray,Chandrasekaran,Wierman,Newton Polytopes and Relative Entropy Optimization.pdf ^NMC20
- 2020,Murray,Chandrasekaran,Wierman,Signomial and polynomial optimization via relative entropy and partial dualization.pdf ^NMC20b
- 2020,Wang,Jaini,Yu, Poupart,A Positivstellensantz for Conditional SAGE Signomials.pdf ^WJYP20
- 2021(thesis), Murray, Applications of convex analysis to signomial and polynomial nonnegativity problems.pdf ^M21
- 2021,Dressler,Murray,Algebraic Perspectives on Signomial Optimization.pdf ^DM21
- 2021,NAUMANN,THEOBALD,SUBLINEAR CIRCUITS FOR POLYHEDRAL SETS.pdf ^NT21 ^DT21

### Slides 
- Thorsten,Relative entropy programming in constrained polynomial ^Thor
- Westerlund,Lundell,Global Optimization of Signomial Programming Problems ^WestLund

---
### Rank/Reason/Summarize

1. [[#^CS14]] This is the paper that started signomial interest  in N&O
	1. Lemma linking AM-GM Expo and relative entropy
	2. Hierarchies (unconstrained (3.3) and constrained (3.6))
	3. Convergence results (unconstrained (Theorem 4.1) and constrained (Theorem 4.2)) Can make pictures based on **extremal exponents**
	4. Dual perspective ?
	5. Collection of why/when sig. opt. is better than POP (discussion section). 
2. ** [[#^DM21]]** Recommended my both M. and L. 
	1. Contains the strongest Positivstellensatz that I know of
	2. Compares to other Pos-stel.
	3. Many nice examples and picture ideas
	4. Tight little proof of NP hardness for sig. opt.
	5. Lots of useful theory.
3. [[#^M21 ]] Hinted by M. and is a thesis on the topic
4. [[#^WJYP20 ]]  Cited a lot by [[#^DM21]]
5. [[#^CS16 ]] Has applications in the title.
6. [[#^MNT20 ]] Murray bias and content could be contained in thesis
7. [[#^NMC20 ]] Murray bias and content could be contained in thesis
8. [[#^NMC20b ]] Murray bias and content could be contained in thesis
10. [[#^KNT20]] Temporal bias
11. [[#^DT21]]
12. [[#^DIW17]]
13. [[#^IW16]]
14. [[#^X14]]
15. [[#^DPZ68 ]] Geometric programming is a specific instance of sig. opt.
16. [[#^C05 ]] Useful only for a singe application example
17. [[#^GIT14 ]] looks tangential
18. [[ ^MHR19]]
19. [[#^MF97]] Useful only for a singe application example
20. [[#^^SW96]]  Useful only for a singe application example


---

#### Applications list:
- From [[#^CS14]]
	- resource allocation in networks [[#^C05 ]] (This looks a bit heavy but there is stuff on page 74)
	- control problems involving chemical processes  [[#^MF97]]  4.2 Alkylation process design
	- spatial frame design [33],  [https://doi.org/10.1016/S0045-7825(96)01097-3]
	- and certain nonlinear flow problems in graphs [21]. [E. L. Peterson, Geometric programming, SIAM Rev., 18 (1976), pp. 1–51]()
- From [[#^DM21]] 
	- Lukas Wachter, Orcun Karaca, Georgios Darivianakis, and Themistoklis Charalambous
	-  chemical reaction networks and biochemical networks [[ ^MHR19]]
	-  And a few more but I think this will have to suffice for now.
	-  




 










--- 
