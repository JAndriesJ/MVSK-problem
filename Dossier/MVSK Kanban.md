---

kanban-plugin: basic

---

## Backlog

- [ ] 8 **Introduction** Literature review on MVSK.
- [ ] 1 **Introduction**: Introduce the MVSK problem
- [ ] 3 **Introduction** Give a mini-survey on the work done on it so-far. (timeline)
- [ ] 2 **Introduction** List all the optimization methods applied.
- [ ] 1 **Introduction** Give the structure of the paper.
- [ ] 1 **Introduction** main contribution of this paper.
- [ ] 1 **Introduction** Tie it all up
- [ ] Write the **Preliminaries** section
- [ ] 3,1 **Preliminaries**: POP results
- [ ] 3,3 **Preliminaries**: Sig-Opt. Results
- [ ] 8,3 **Preliminaries**: Other optimization
- [ ] 1,1 Round off the **Modeling** section
- [ ] Collect **numerical results** and experiments in a section
- [ ] **Numerical:** Compare POP bounds for different levels.
- [ ] **Numerical:** Compare bounds for different opt. techniques.
- [ ] **Numerical:** Compare bounds for different levels of sig. opt.
- [ ] **Numerical:** Compare bounds with and without card restrictions.
- [ ] **Numerical:** Visualize the data.
- [ ] Write a section on **Extensions/Discussion** and further ideas
- [ ] Question: do higher order moments get drowned out?
- [ ] Formalize the effect of cardinality constraints on POP, Sig.Opt., and other methods.
- [ ] Is there a correlation between higher order moment behavior and price movement?
- [ ] Is weak SDP optimization viable here?
- [ ] When does past data cease to be relevant? When does the echo fade to static?
- [ ] Does sparsity play a role?
- [ ] Make a dash board that tracks the moments of stocks


## Priority

- [ ] 5,8 Make an easy to use **demo for computing** numerical results. **(gold)**
- [ ] ?,5 **Code**: Install NAG library
- [ ] ?,?**Code**: minimal working example of multi-start
- [ ] ?,? **Code:** Minimal example of particle swarm
- [ ] 3,3 **Code**: Make a readme.txt


## In Progress

- [ ] 3,3**Code**: get the Sig-Opt. to work inside of Julia. Focus just on the first two moments.


## Done

**Complete**
- [x] **Code:** Data_moments module
- [x] 2,5 **Code**: Find alternative POP code and make a minimal example.
- [x] 1,1 **Code**: Modularize the convex Quadratic optimization. ()~2 hours)


## Things to remember

- [ ] I.N.V.E.S.T. (Independent, Negotiable, Valuable, Estimable, Small, Timely)
- [ ] Ranked 1,2,3,5,8,13,21
- [ ] Goal are narrative.
- [ ] Demo or Die.
- [ ] Immediate error correction
- [ ] Do nothing that can be delegated.
- [ ] Assume a state of liminality.


## Stopped

- [ ] 3,2 **Preliminaries** Cardinality constraints
- [ ] 2,2**Code**:  Fix the POP code




%% kanban:settings
```
{"kanban-plugin":"basic"}
```
%%