# CLL_IRL
Stable Inverse Reinforcement Learning: Policies from Control Lyapunov Landscapes
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
## Usage

There are two notebooks provided for single shape evaluation or multi shape evaluation. They should run out of the box. It is assumed that notebooks are executed while being in the root-path of this repository (otherwise the addpath command will not work)

---
## Structure

### functions
Includes the algorithm itself (func_switching_stochastic_gradsquare.mlx), a plotting and simulation utility (func_simulate_clfsos.mlx) and a dataloading utility (plot_shape.m)

### LASA
The LASA dataset.

### misc
Some helpful stuff

### notebooks
Two example notebooks how to use the algorithm
