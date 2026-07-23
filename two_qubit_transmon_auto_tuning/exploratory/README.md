# Exploratory / prototyping scripts

These are the early, standalone scripts written while building up to the main
simulation pipeline. They are kept for provenance — they show how the model was
developed from a single qubit to two coupled transmons and how the optimization
approach was validated — but they are **not** part of the maintained package and
may use earlier conventions or be incomplete.

The maintained, tested pipeline lives in the top-level package
(`simulation.py`, `optimization.py`).

| File | What it explores |
|------|------------------|
| `single_transmon.py` | Time evolution of one transmon; sanity-checks the drift Hamiltonian and QuTiP set-up. |
| `two_transmons.py` | Two coupled transmons evolved from several initial states, visualized on the Bloch sphere. |
| `qpt_2level_prototype.py` | First quantum-process-tomography pipeline on a 2-level model (before leakage was added), including power-Rabi and an X-gate fidelity check. |
| `bayesian_from_scratch.py` | A from-scratch Gaussian-process Bayesian optimizer (scikit-learn GP + probability-of-improvement acquisition) used to understand the method before adopting the `bayesian-optimization` library. |

Each script is self-contained; run it directly with the project environment
active, e.g. `python two_qubit_transmon_auto_tuning/exploratory/single_transmon.py`.
