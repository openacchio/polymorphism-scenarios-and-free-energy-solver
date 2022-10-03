# Polymorphism-scenarios-and-free-energy-solver
This repository contains Matlab and Python code to reproduce the results of the paper "Multifarious selection on the evolution of visual signals", <to be added here: ref of paper; link: TBD>. 

The code tracks the evolution of the distribution of phenotypes in a population with a constant number of individuals in different scenarios. Depending on the different forces at play operate synergistically or in opposition, the distribution of phenotypes evolves towards all types of phenotypic variation, including monomorphism, continuous variation, and discrete polymorphism.  

The Matlab function `create_landscape_for_scenarios_ABC.m` generates the three components that make the fitness landscape, namely F_cost, F_nat and F_sex, and exports them to build the total fitness landscape F_total = F_cost + F_nat + F_sex.

The Python function `Pop_Evo_multifarious_evolution_free_energy_Fokker_Planck_equation.py` models the temporal evolution of the distribution of phenotypes given a random component (the diffusion term, D) and the fitness landscape F_total. This is done using a Fokker-Planck equation in its equivalent formulation using free energy equivalent (see text for details). The code finds the distribution of phenotypes that minimizes free energy, which corresponds to the distribution that maximizes the fitness on the fitness landscape F_total given the processes of diffusion driven by D.  

The evolutionary trajectories of the distribution of phenotypes for all scenarios are provided by videos `movie_scenarioA_compressed.mp4`, `movie_scenarioB1_compressed.mp4`, `movie_scenarioC_case1_compressed.mp4`, and `movie_scenarioC_case2_compressed.mp4`.
