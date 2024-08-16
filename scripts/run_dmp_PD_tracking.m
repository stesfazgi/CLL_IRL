clear;
clc;
yalmip('clear');
% addpath('..\functions');
addpath('..\LASA\DataSet\');
% addpath('..\functions\CLFDM_lib\');
addpath(genpath('..\functions'))

%%%%%%%%%% Evaluation specific parameters
% Controller
P_gain = 200;
D_gain = 1;

% Simulation
n_sim_steps = 2000;

% Input directory (DMP run)
input_dir = '..\eval_output\eval_dmp';

% Output directory
output_algoname = 'dmp_PD';
output_dir = '..\eval_output\eval_dmp_PD';

if ~exist(output_dir, "dir")
    mkdir(output_dir);
    fprintf("Making new directory: %s \n", output_dir);
else
    warning("Run already exists. Continuing...");
end

eval_shape_ids = [1 3 19 23 26 30];

for esid=1:length(eval_shape_ids)
    
    shape_id = eval_shape_ids(esid);

    [~, ~, shapename, dt] = plot_shape(shape_id, 7, 1000, false, [0 0]);
    dmp_plot_data_path = fullfile(input_dir, strcat('plot_data_dmp_', shapename, '.mat'));
    load(dmp_plot_data_path);


    %%%%%%%%% Load unstable system dynamics
    x1 = sdpvar(1, 1);
    x2 = sdpvar(1, 1);
    % State vector
    x = [x1; x2];

    savepath = fullfile("..\eval_sys_unstable", strcat(string(shape_id), ".mat"));
    if ~exist(savepath, "file")
        error("System not found!")
    else
        [f, f_coeffs] = func_restore_poly_sys(x, savepath);
        fprintf("Using old system.\n");
    end

    % Input matrix of system
    g = eye(2);

    algo_params.P_gain = P_gain;
    algo_params.D_gain = D_gain;

    % Load DMP model
    in = 1:algo_params.Vxf.d;
    out = algo_params.Vxf.d+1:2*algo_params.Vxf.d;
    fn_handle_GMR = @(x) GMR(algo_params.Priors_EM, algo_params.Mu_EM, algo_params.Sigma_EM, x, in,out);
    % Get safe DS handle
    fn_handle = @(x) DS_stabilizer(x,fn_handle_GMR,algo_params.Vxf,algo_params.rho0,algo_params.kappa0);

    reproductions = {};

    for t=1:data_params.n_demos
        curr_repro = zeros(2, n_sim_steps + 1);
        x_curr = target_data{t}(:, 1);
        curr_repro(:, 1) = x_curr;

        for i=1:n_sim_steps+1
            
            x_des = x_curr + dt*fn_handle(x_curr);
            dx_des = fn_handle(x_curr);
            dx_curr = replace(f, x, x_curr);
        
            e_x = (x_des - x_curr);
            e_dx = (dx_des - dx_curr);
        
            u = P_gain*e_x + D_gain*e_dx;
        
            dx = dx_curr + g*u;
        
            x_curr = x_curr + dx*dt;
            curr_repro(:, i) = x_curr;

        end

        reproductions{end+1} = curr_repro;

    end

    algoname = output_algoname;

    output_file = fullfile(output_dir, strcat('plot_data_', algoname, '_', shapename, '.mat'));
    if exist(output_file, "file")
        warning("Shape %s already exist. Skipping...", shapename);
        continue;
    end

    save(output_file, 'shapename', 'algoname', 'algo_params', 'data_params', ...
        'v_surf', 'x1_dot_surf', 'x2_dot_surf', 'reproductions', ...
        'x1_range', 'x2_range', 't_algo', 'target_data');
    
end


beep;

