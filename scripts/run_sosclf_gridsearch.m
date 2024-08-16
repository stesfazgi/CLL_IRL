%% Script to perform a gridsearch over a number of parameters
clear;
clc;
yalmip('clear');
addpath('..\functions');
addpath('..\LASA\DataSet\');

gs_rootfolder = '..\gridsearch_output\gridsearch13_1000_highDegB';

% Simulation/Plotting settings
plot_padding = 0.5;
n_interpolation = 100;
n_sim_steps = 2000;

% Number of demonstrations to use
n_demos = 7;

% Values for degree of V
gs_degV = [8];
% Values for degree of B
gs_degB = [16];
% Number of samples per demonstration
gs_n_samples = [1000];

% Possible eigendynamics of the system
% 0=No eigendynamics (f=0)
% 1=Instable eigendynamics
gs_eigendyns = [1];

% Weighting factor for regularization term (during V and B optimization
% respectively)
gs_regV = [1e-6, 1e-4, 1e-3];
gs_regB = [1e-6, 1e-4, 1e-3];

% Number of iterations
gs_n_iter = [5];

% IDs of shapes to evaluate
%   ID      NAME                ID      Name            ID      Name
%   1       Angle               11      Leaf_2          21      Spoon
%   2       BendedLine          12      Line            22      Sshape
%   3       CShape              13      LShape          23      Trapezoid
%   4       DoubleBendedLine    14      NShape          24      Worm
%   5       GShape              15      PShape          25      WShape
%   6       heee                16      RShape          26      Zshape
%   7       JShape              17      Saeghe          27      Multi_Models_1
%   8       JShape_2            18      Sharpc          28      Multi_Models_2
%   9       Khamesh             19      Sine            29      Multi_Models_3
%  10       Leaf_1              20      Snake           30      Multi_Models_4
gs_shapeids = [1 3 19 23 26 30];
% gs_shapeids = [1 26];
% gs_shapeids = linspace(1, 30, 30);

% Create search grid
[DV, DB, NS, ED, RV, RB, NI] = ndgrid(gs_degV, gs_degB, gs_n_samples, ...
    gs_eigendyns, gs_regV, gs_regB, gs_n_iter);

% Create rootdirectory and save configuration
if ~exist(gs_rootfolder, "dir")
    mkdir(gs_rootfolder);
    n_runs = length(DV(:));
    save(fullfile(gs_rootfolder, "root_config.mat"));
else
    warning("Resuming experiment...");
%     error("Gridsearch directory %s already exists. Preventing overwrite.", gs_rootfolder);
end

all_configs = {};

for i=1:length(DV(:))
    cfg.degV = DV(i);
    cfg.degB = DB(i);
    cfg.n_samples = NS(i);
    cfg.eigendyns = ED(i);
    cfg.regV = RV(i);
    cfg.regB = RB(i);
    cfg.n_iter = NI(i);
    cfg.shapeids = gs_shapeids;
    cfg.plot_padding = plot_padding;
    cfg.n_interpolation = n_interpolation;
    cfg.n_sim_steps = n_sim_steps;
    cfg.n_demos = n_demos;
    % Create unique label for run
    run_name = sprintf("eval_%i", i);
    cfg.root = fullfile(gs_rootfolder, run_name);
    all_configs{end+1} = cfg;
end

fprintf('%i possible combination(s) found. \n', length(all_configs));

t_eval_start = tic;
%% Run gridsearch

% c1 = parcluster;

% myjobs = {};
% for i=1:length(all_configs)
%     j = batch(@func_run_sosclf_gs_reg_from_config, 1, {all_configs{i}}, 'Pool', 1);
% end
% wait(j);
for i=1:length(all_configs)
    func_run_sosclf_gs_reg_from_config(all_configs{i});
end
t_eval_end = toc(t_eval_start);
fprintf("Took %0.3f seconds\n", t_eval_end);
beep;