
# Two-Qubit-Transmon-Auto-Tuning
Codespace to simulate superconducting qubit gates as well as machine learning algorithms to optimize them.

qpt_simulation.py is a script that simulates the 2-qubit setup, returning the chi-matrix for a quantum process defined by some defined drive pulses.

fidelity_function.py nests the above script into functional form, interfacing it with the Bayesian optimization code. Upon defining the pulse parameters, the 
fidelity with respect to some chosen gate is returned. It was created by transforming qpt_simulation.py directly, and should be refined to less complexity
and higher readability.

bayesian_stage1.py runs the Bayesian optimization procedure, for whichever gate is defined as the objective.

bayesian_stage0.py was an attempt to understand the inner workings of the Bayesian process.

