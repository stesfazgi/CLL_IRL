%% This script runs the ground truth evaluation and save results as plotting data

clc;
clear;
close all;
yalmip('clear');
cd('..');
addpath(genpath('functions\'));

% ---------------------- ALGORITHM PARAMS ----------------------
% Degree of V
degV = 8;
% Degree of B
degB = 2;
% Regularization factor for fitting V
regV = 1e-4;
% Regularization factor for fitting B
regB = 1e-4;
% Number of iterations for fitting V
n_iter = 3;
% Coefficient tolerance. Parameters with magnitude smaller than "tol" are
% set to 0
tol = 1e-6;

% ---------------------- PLOTTING PARAMS ----------------------
% Mesh size per axis for evaluation
n_mesh = 100;
% Creates a figure showning the results in matlab (for debug)
create_figures = true;
% Path under which to store the output data
output_path = "eval_latex/data_gt";






if ~exist(output_path, "dir")
    mkdir(output_path);
else
    error("Output path %s already exists. Preventing overwrite.", output_path);
end

x_range = linspace(-1, 1, n_mesh);
[Xte1, Xte2] = meshgrid(x_range, x_range);
Xte = [Xte1(:), Xte2(:)];


% Seed random number generator
rng(1);

% Define system states and ground truth lyapunov (V*)
x = sdpvar(2, 1);
V_star = 0.5*x(2)^2 + 0.1*x(1)^2 + x(1)*x(2)^2 + x(1)^4 + x(2)^4;
dVdx_star = jacobian(V_star, x);

% Verify that V* is SOS
res = solvesos(sos(V_star), [], sdpsettings('verbose', 0));
if res.problem == 0
    fprintf("V* is SOS and fulfills first lyapunov constraint. \n");
else
    error("V is not SOS")
end

% Verify that V* fulfills proxy constraint
res2 = solvesos(sos(dVdx_star*dVdx_star' - 1e-3*sum(x.^2)), [], sdpsettings('verbose', 0));
if res2.problem == 0
    fprintf("dVdx* fulfills proxy constraint. \n");
else
    error("GS is not SOS")
end

% Define g
g = eye(2);

% Define f
degF = 3;
f = sdpvar(2, 1);
[f(1), f1_coeffs, f1z] = polynomial(x, degF, 1);
[f(2), f2_coeffs, f2z] = polynomial(x, degF, 1);

% Define Q (state costs)
[Q, Q_coeffs, Qm] = polynomial(x, 4, 2);

% Define the HJB that has to be solved
HJB = dVdx_star*f - 0.25*dVdx_star*g*g'*dVdx_star' + Q;

% Solve for coefficients of f
[HJB_coeffs, ~] = coefficients(HJB, x);
res = optimize([HJB_coeffs == 0, sos(Q - 1e-3*sum(x.^2))], [], sdpsettings('verbose', 0), [f1_coeffs; f2_coeffs; Q_coeffs]);
if res.problem == 0
    fprintf("Successfully solved for f and Q.\n");
else
    fprintf("No solution found. Problem %i\n", res.problem);
end

% Get coefficients for f and remove numerically insignificant terms
f1_coeffs_clean = (abs(value(f1_coeffs)) > tol).*value(f1_coeffs);
f2_coeffs_clean = (abs(value(f2_coeffs)) > tol).*value(f2_coeffs);

% Fitted polynomial system
f_star = replace(f, [f1_coeffs; f2_coeffs], [f1_coeffs_clean; f2_coeffs_clean]);

% Closed loop system under V*
f_cl_star = f_star - 0.5*g*g'*dVdx_star';

% Verify that stability of f_cl_star is certified by V_star
res = solvesos(sos(-dVdx_star*f_cl_star), [], sdpsettings('verbose', 0), []);
if res.problem == 0
    fprintf("f* and V* fulfill second lyapunov constraint.\n");
else
    error("Second lyapunov constraint violated for f* and V*!")
end

% Fitted symbolic lyapunov function
% Reconstruct polynomial
[Vhat, Vhat_coeffs, Vhat_m] = polynomial(x, degV, 2);
dVdx_hat = jacobian(Vhat, x);

% Fitted symbolic beta function
[B1hat, B1hat_coeffs, ~] = polynomial(x, degB, 0);
[B2hat, B2hat_coeffs, ~] = polynomial(x, degB, 0);

% Expression of closed loop dynamics depending on x and Vhat_coeffs
f_cl = f_star - diag([B1hat, B2hat])*g*g'*dVdx_hat';

% Initial guess for lyapunov function
% Note: only for plotting, only init_scale is passed to the algorithm which
% automatically create the quadratic form.
init_scale = 0.1;
V_init = init_scale*(x(1)^2 + x(2)^2);
dVdx_init = jacobian(V_init, x);
% Expression of closed loop dynamics under V_init
f_cl_init = f_star - 0.5*g*g'*dVdx_init';

% Create options struct
opt.degV = degV;
opt.degB = degB;
opt.gammaV = regV;
opt.gammaB = regB;
opt.n_iter = n_iter;
opt.init_scale = init_scale;


% ---------------------- EXPERIMENT 1 ----------------------
% 7x7 gridpoints in interval [-1, 0] for each state
x1_bounds = [-1 0];
x2_bounds = [-1 0];
nx1_samples = 7;
nx2_samples = 7;

sigma = 0.001*(sum(x.^2));

data_pos = [];
data_vel = [];

x1_samp = linspace(x1_bounds(1), x1_bounds(2), nx1_samples);
x2_samp = linspace(x2_bounds(1), x2_bounds(2), nx2_samples);
for i=1:nx1_samples
    for j=1:nx2_samples
        x_curr = [x1_samp(i); x2_samp(j)];
        % Add noise to target dynamics
        x_dot = replace(f_cl_star, x, x_curr) + randn(2, 1)*sqrt(replace(sigma, x, x_curr));
        data_pos = [data_pos x_curr];
        data_vel = [data_vel x_dot];
    end
end
% Fit system
[Vc_hat_exp1, B1c_hat_exp1, B2c_hat_exp1] = func_sosclf_gs_reg(x, f_star, g, opt, data_pos, data_vel);
% Clean coefficients
Vc_hat_exp1 = (abs(Vc_hat_exp1) > tol).*Vc_hat_exp1;
% Store data in experiment struct
exp1.target_data = {data_pos};



% ---------------------- EXPERIMENT 2 ----------------------
% 7x14 gridpoints in interval [-1, 0] for x1 and [-1, 1] for x2
x1_bounds = [-1 0];
x2_bounds = [-1 1];
nx1_samples = 7;
nx2_samples = 14;

sigma = 0.001*(sum(x.^2));

data_pos = [];
data_vel = [];

x1_samp = linspace(x1_bounds(1), x1_bounds(2), nx1_samples);
x2_samp = linspace(x2_bounds(1), x2_bounds(2), nx2_samples);
for i=1:nx1_samples
    for j=1:nx2_samples
        x_curr = [x1_samp(i); x2_samp(j)];
        % Add noise to target dynamics
        x_dot = replace(f_cl_star, x, x_curr) + randn(2, 1)*sqrt(replace(sigma, x, x_curr));
        data_pos = [data_pos x_curr];
        data_vel = [data_vel x_dot];
    end
end
% Fit system
[Vc_hat_exp2, B1c_hat_exp2, B2c_hat_exp2] = func_sosclf_gs_reg(x, f_star, g, opt, data_pos, data_vel);
% Clean coefficients
Vc_hat_exp2 = (abs(Vc_hat_exp2) > tol).*Vc_hat_exp2;
% Store data in experiment struct
exp2.target_data = {data_pos};



% ---------------------- EXPERIMENT 3 ----------------------
% 14x14 gridpoints in interval [-1, 1] for each state
x1_bounds = [-1 1];
x2_bounds = [-1 1];
nx1_samples = 14;
nx2_samples = 14;

sigma = 0.001*(sum(x.^2));

data_pos = [];
data_vel = [];

x1_samp = linspace(x1_bounds(1), x1_bounds(2), nx1_samples);
x2_samp = linspace(x2_bounds(1), x2_bounds(2), nx2_samples);
for i=1:nx1_samples
    for j=1:nx2_samples
        x_curr = [x1_samp(i); x2_samp(j)];
        % Add noise to target dynamics
        x_dot = replace(f_cl_star, x, x_curr) + randn(2, 1)*sqrt(replace(sigma, x, x_curr));
        data_pos = [data_pos x_curr];
        data_vel = [data_vel x_dot];
    end
end
% Fit system
[Vc_hat_exp3, B1c_hat_exp3, B2c_hat_exp3] = func_sosclf_gs_reg(x, f_star, g, opt, data_pos, data_vel);
% Clean coefficients
Vc_hat_exp3 = (abs(Vc_hat_exp3) > tol).*Vc_hat_exp3;
% Store data in experiment struct
exp3.target_data = {data_pos};



% ----------------------- PLOTTING -----------------------
% Get symbolic expressions for estimated V in each experiment
V_hat_exp1 = replace(Vhat, Vhat_coeffs, Vc_hat_exp1);
V_hat_exp2 = replace(Vhat, Vhat_coeffs, Vc_hat_exp2);
V_hat_exp3 = replace(Vhat, Vhat_coeffs, Vc_hat_exp3);
% Get symbolic expressions for closed loop systems
f_cl_exp1 = replace(f_cl, [Vhat_coeffs; B1hat_coeffs; B2hat_coeffs], [Vc_hat_exp1; B1c_hat_exp1; B2c_hat_exp1]);
f_cl_exp2 = replace(f_cl, [Vhat_coeffs; B1hat_coeffs; B2hat_coeffs], [Vc_hat_exp2; B1c_hat_exp2; B2c_hat_exp2]);
f_cl_exp3 = replace(f_cl, [Vhat_coeffs; B1hat_coeffs; B2hat_coeffs], [Vc_hat_exp3; B1c_hat_exp3; B2c_hat_exp3]);

% Transform into function handles for faster evaluation over grid
h_V_hat_exp1 = func_sdpvar2fnhandle(V_hat_exp1, x);
h_V_hat_exp2 = func_sdpvar2fnhandle(V_hat_exp2, x);
h_V_hat_exp3 = func_sdpvar2fnhandle(V_hat_exp3, x);
% Also transform initial guess and ground truth
h_V_init = func_sdpvar2fnhandle(V_init, x);
h_V_star = func_sdpvar2fnhandle(V_star, x);

h_f_cl_exp1 = func_sdpvar2fnhandle(f_cl_exp1, x);
h_f_cl_exp2 = func_sdpvar2fnhandle(f_cl_exp2, x);
h_f_cl_exp3 = func_sdpvar2fnhandle(f_cl_exp3, x);

h_f_ol = func_sdpvar2fnhandle(f_star, x);
h_f_cl_star = func_sdpvar2fnhandle(f_cl_star, x);
h_f_cl_init = func_sdpvar2fnhandle(f_cl_init, x);

% Initialize meshes to hold lyapunov evaluations on grid
[v_surf_init, v_surf_star, v_surf_exp1, v_surf_exp2, v_surf_exp3] = deal(zeros(n_mesh, n_mesh));

% Initialize tensors to hold dynamics evaluated on grid
[dx_ol, dx_cl_init, dx_cl_star, dx_cl_exp1, dx_cl_exp2, dx_cl_exp3] = deal(zeros(n_mesh, n_mesh, 2));

% Create plotting data
for i=1:n_mesh
    for j=1:n_mesh
        x_curr = [x_range(i); x_range(j)];
        
        % Evaluate dynamics
        dx_ol(j, i, :) = h_f_ol(x_curr);
        dx_cl_init(j, i, :) = h_f_cl_init(x_curr);
        dx_cl_star(j, i, :) = h_f_cl_star(x_curr);
        dx_cl_exp1(j, i, :) = h_f_cl_exp1(x_curr);
        dx_cl_exp2(j, i, :) = h_f_cl_exp2(x_curr);
        dx_cl_exp3(j, i, :) = h_f_cl_exp3(x_curr);

        % Evaluate lyapunov/value functions
        v_surf_init(j, i) = h_V_init(x_curr);
        v_surf_star(j, i) = h_V_star(x_curr);
        v_surf_exp1(j, i) = h_V_hat_exp1(x_curr);
        v_surf_exp2(j, i) = h_V_hat_exp2(x_curr);
        v_surf_exp3(j, i) = h_V_hat_exp3(x_curr);
    end
end

% Put plotting data for each experiment in separate structs
algo_params.degV = degV;
algo_params.degB = degB;
algo_params.n_iter = n_iter;
algoname = 'sosclf';

exp1.algo_params = algo_params;
exp1.algoname = algoname;
exp1.x1_range = x_range;
exp1.x2_range = x_range;
exp1.v_surf = v_surf_exp1;
exp1.v_surf_star = v_surf_star;
exp1.x1_dot_surf = dx_cl_exp1(:, :, 1);
exp1.x2_dot_surf = dx_cl_exp1(:, :, 2);

exp2.algo_params = algo_params;
exp2.algoname = algoname;
exp2.x1_range = x_range;
exp2.x2_range = x_range;
exp2.v_surf = v_surf_exp2;
exp2.v_surf_star = v_surf_star;
exp2.x1_dot_surf = dx_cl_exp2(:, :, 1);
exp2.x2_dot_surf = dx_cl_exp2(:, :, 2);

exp3.algo_params = algo_params;
exp3.algoname = algoname;
exp3.x1_range = x_range;
exp3.x2_range = x_range;
exp3.v_surf = v_surf_exp3;
exp3.v_surf_star = v_surf_star;
exp3.x1_dot_surf = dx_cl_exp3(:, :, 1);
exp3.x2_dot_surf = dx_cl_exp3(:, :, 2);

exp_init.algo_params = algo_params;
exp_init.algoname = algoname;
exp_init.x1_range = x_range;
exp_init.x2_range = x_range;
exp_init.v_surf = v_surf_init;
exp_init.v_surf_star = v_surf_star;
exp_init.x1_dot_surf = dx_cl_init(:, :, 1);
exp_init.x2_dot_surf = dx_cl_init(:, :, 2);

exp_star.algo_params = algo_params;
exp_star.algoname = algoname;
exp_star.x1_range = x_range;
exp_star.x2_range = x_range;
exp_star.v_surf = v_surf_star;
exp_star.v_surf_star = v_surf_star;
exp_star.x1_dot_surf = dx_cl_star(:, :, 1);
exp_star.x2_dot_surf = dx_cl_star(:, :, 2);

exp_ol.algo_params = algo_params;
exp_ol.algoname = algoname;
exp_ol.x1_range = x_range;
exp_ol.x2_range = x_range;
exp_ol.x1_dot_surf = dx_ol(:, :, 1);
exp_ol.x2_dot_surf = dx_ol(:, :, 2);

% Export plotting data
save(fullfile(output_path, "plot_data_gt_eval_init.mat"), '-struct', "exp_init");
save(fullfile(output_path, "plot_data_gt_eval_star.mat"), '-struct', "exp_star");
save(fullfile(output_path, "plot_data_gt_eval_ol.mat"), '-struct', "exp_ol");
save(fullfile(output_path, "plot_data_gt_eval1.mat"), '-struct', "exp1");
save(fullfile(output_path, "plot_data_gt_eval2.mat"), '-struct', "exp2");
save(fullfile(output_path, "plot_data_gt_eval3.mat"), '-struct', "exp3");

if create_figures
    clfsos_levels = logspace(-2, log10(max(v_surf_star(:))), 10);

    figure('Position', [0, 0, 1500, 600]);
    hold on;
    % STAR
    subplot(2, 5, 1);
    l = streamslice(x_range, x_range, dx_cl_star(:, :, 1), dx_cl_star(:, :, 2));
    hold on;
    set(l, 'LineWidth', 0.8);
    set(l, 'Color', [.0 .0 .0]);
    set(l, 'HandleVisibility', 'off');
    xlim([-1, 1]);
    ylim([-1, 1]);
    title("Ground Truth");
    subplot(2, 5, 6);
    contour(x_range, x_range, v_surf_star, clfsos_levels, 'b--', 'LineWidth', 1.5, 'DisplayName', 'V star');
    xlim([-1, 1]);
    ylim([-1, 1]);
    legend;

    % INIT
    subplot(2, 5, 2);
    l = streamslice(x_range, x_range, dx_cl_init(:, :, 1), dx_cl_init(:, :, 2));
    hold on;
    set(l, 'LineWidth', 0.8);
    set(l, 'Color', [.0 .0 .0]);
    set(l, 'HandleVisibility', 'off');
    xlim([-1, 1]);
    ylim([-1, 1]);
    title("Initial guess");
    subplot(2, 5, 7);
    contour(x_range, x_range, v_surf_init, clfsos_levels, 'g', 'LineWidth', 1.5, 'DisplayName', 'V_{init}');
    hold on;
    contour(x_range, x_range, v_surf_star, clfsos_levels, 'b--', 'LineWidth', 1.5, 'DisplayName', 'V star');
    xlim([-1, 1]);
    ylim([-1, 1]);
    legend;

    % EXP1
    subplot(2, 5, 3);
    l = streamslice(x_range, x_range, dx_cl_exp1(:, :, 1), dx_cl_exp1(:, :, 2));
    hold on;
    set(l, 'LineWidth', 0.8);
    set(l, 'Color', [.0 .0 .0]);
    set(l, 'HandleVisibility', 'off');
    xlim([-1, 1]);
    ylim([-1, 1]);
    title("Experiment 1");
    subplot(2, 5, 8);
    contour(x_range, x_range, v_surf_exp1, clfsos_levels, 'g', 'LineWidth', 1.5, 'DisplayName', 'V_1');
    hold on;
    contour(x_range, x_range, v_surf_star, clfsos_levels, 'b--', 'LineWidth', 1.5, 'DisplayName', 'V star');
    plot(exp1.target_data{1}(1, :), exp1.target_data{1}(2, :), 'r+', 'LineWidth', 1.5, 'DisplayName', 'Data');
    xlim([-1, 1]);
    ylim([-1, 1]);
    legend;

    % EXP2
    subplot(2, 5, 4);
    l = streamslice(x_range, x_range, dx_cl_exp2(:, :, 1), dx_cl_exp2(:, :, 2));
    hold on;
    set(l, 'LineWidth', 0.8);
    set(l, 'Color', [.0 .0 .0]);
    set(l, 'HandleVisibility', 'off');
    xlim([-1, 1]);
    ylim([-1, 1]);
    title("Experiment 2");
    subplot(2, 5, 9);
    contour(x_range, x_range, v_surf_exp2, clfsos_levels, 'g', 'LineWidth', 1.5, 'DisplayName', 'V_2');
    hold on;
    contour(x_range, x_range, v_surf_star, clfsos_levels, 'b--', 'LineWidth', 1.5, 'DisplayName', 'V star');
    plot(exp2.target_data{1}(1, :), exp2.target_data{1}(2, :), 'r+', 'LineWidth', 1.5, 'DisplayName', 'Data');
    xlim([-1, 1]);
    ylim([-1, 1]);
    legend;

    % EXP3
    subplot(2, 5, 5);
    l = streamslice(x_range, x_range, dx_cl_exp3(:, :, 1), dx_cl_exp3(:, :, 2));
    hold on;
    set(l, 'LineWidth', 0.8);
    set(l, 'Color', [.0 .0 .0]);
    set(l, 'HandleVisibility', 'off');
    xlim([-1, 1]);
    ylim([-1, 1]);
    title("Experiment 3");
    subplot(2, 5, 10);
    contour(x_range, x_range, v_surf_exp3, clfsos_levels, 'g', 'LineWidth', 1.5, 'DisplayName', 'V_3');
    hold on;
    contour(x_range, x_range, v_surf_star, clfsos_levels, 'b--', 'LineWidth', 1.5, 'DisplayName', 'V star');
    plot(exp3.target_data{1}(1, :), exp3.target_data{1}(2, :), 'r+', 'LineWidth', 1.5, 'DisplayName', 'Data');
    xlim([-1, 1]);
    ylim([-1, 1]);
    legend;

end
