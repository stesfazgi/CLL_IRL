# Stable Inverse Reinforcement Learning: Policies from Control Lyapunov Landscapes

This repository contains an implementation of the paper [1]:

[1] Tesfazgi, S., Sprandl, L., Lederer, A., & Hirche, S. "Stable Inverse Reinforcement Learning: Policies from Control Lyapunov Landscapes" in IEEE Open Journal of Control Systems, 2024 [full paper](https://ieeexplore.ieee.org/document/10643266)

---

## Requirements
- Matlab (tested on R2021a)
- Mosek (Follow the instructions at https://docs.mosek.com/10.0/toolbox/install-interface.html)
- Yalmip (Follow the instructions at https://yalmip.github.io/tutorial/installation/)
- LASA dataset (create a folder called LASA and include contents from https://bitbucket.org/khansari/lasahandwritingdataset/src/master/)
---
## Installation
1. Install SeDuMi and then Yalmip by following their respective installation instructions.
2. (Optionally) Include the path to SeDuMi and Yalmip into your matlab startup file. This way you dont need to rerun the installation commands at each Matlab startup manually:

    2.1. Open your Matlab startup.m file by running: ```edit(fullfile(userpath,'startup.m'))```

    2.2. Append the contents of **misc/startup.m** to your startup file and adjust the Yalmip and sedumi path to your installation path.

---
## Structure

### functions
Includes the algorithm itself (func_switching_stochastic_gradsquare.mlx), a plotting and simulation utility (func_simulate_clfsos.mlx) and a dataloading utility (plot_shape.m)

### LASA
The LASA dataset.

### misc
Some helpful stuff

### scripts
main code including example uses of the algorithm

---
## References
If you found this software useful for your research, consider citing us.
```
@ARTICLE{Tesfazgi24,
  author={Tesfazgi, Samuel and Sprandl, Leonhard and Lederer, Armin and Hirche, Sandra},
  journal={IEEE Open Journal of Control Systems}, 
  title={Stable Inverse Reinforcement Learning: Policies from Control Lyapunov Landscapes}, 
  year={2024},
  volume={},
  number={},
  pages={1-17},
  doi={10.1109/OJCSYS.2024.3447464}}
```
