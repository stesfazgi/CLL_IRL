clear;
clc;

yalmip('clear');
addpath(genpath('..\functions'));
addpath('..\LASA\DataSet\');

%%%%%%%%%% Evaluation specific parameters
% Dataset
n_samples = 1000;
n_demos = 7;

% Simulation/Plotting
plot_padding = 0.5;
n_interpolation = 100;
n_sim_steps = 2000;

% Output directory
output_dir = '..\eval_output\eval_dmp\';

if exist(output_dir, "dir")
    error("DMP output dir already exists. Preventing overwriting: %s", output_dir);
else
    mkdir(output_dir);
end

% Shape-IDs to consider
shape_ids = [1 3 19 23 26 30];
shape_seed = [9 1 15 13 13 15];

%%%%%%%%%% Algorithm specific parameters
% Dimension of state (default)
Vxf0.d = 2;
% Objective weighting
Vxf0.w = 1e-4;
% Number of WSAQF components
Vxf0.L = 3;

%%%%%%%%%% Algo Init
rand_init = true;

options.tol_mat_bias = 10^-15;       
options.display = 1;   
options.tol_stopping=10^-10;
options.max_iter = 1000;
options.optimizePriors = true;
options.upperBoundEigenValue = true;

% Parameters for bounding function that determines when a corrective
% u* has to be applied to force system along decreasing V direction
% Bounding to enforce decrease rate
rho0 = 1e-3;
kappa0 = 1;
% Number of GMM components
n_gmm_components = 8;

    
for curr_id=1:length(shape_ids)

    seed = shape_seed(curr_id);
     % Seed random number generator for initial guesses
    rng(seed);

    shape_id = shape_ids(curr_id);

    %%%%%%%%% Perpare dataset
    % Load dataset
    [data_pos, data_vel, shapename, dt] = plot_shape(shape_id, n_demos, n_samples, false, [0 0]);
    % Normalize dataset
    scale_factor = 1/max(abs(data_pos(:)));
    data_pos = data_pos*scale_factor;
    data_vel = data_vel*scale_factor;

    data_start_idx = linspace(n_samples, n_samples*n_demos, n_demos) - (n_samples - 1);

    clfdm_data = [data_pos; data_vel];

    %%%%%%%%% Plot settings
    axis_bounds = zeros(2, 2);
    axis_bounds(1, :) = [min(data_pos(1, :))-plot_padding max(data_pos(1, :))+plot_padding];
    axis_bounds(2, :) = [min(data_pos(2, :))-plot_padding max(data_pos(2, :))+plot_padding];

    %%%%%%%%% Initial guess
    if rand_init
        lengthScale = sqrt(var(clfdm_data(1:Vxf0.d,:)'));
        lengthScaleMatrix = sqrtm(cov(clfdm_data(1:Vxf0.d,:)'));
        lengthScale = lengthScale(:);
        Vxf0.Priors = rand(Vxf0.L+1,1);
        for l=1:Vxf0.L+1
            tmpMat = randn(Vxf0.d,Vxf0.d);
            Vxf0.Mu(:,l) = randn(Vxf0.d,1).*lengthScale;
            Vxf0.P(:,:,l) = lengthScaleMatrix*(tmpMat*tmpMat')*lengthScaleMatrix;
        end
    else
        Vxf0.Priors = ones(Vxf0.L+1,1);
        Vxf0.Priors = Vxf0.Priors/sum(Vxf0.Priors);
        Vxf0.Mu = zeros(Vxf0.d,Vxf0.L+1);
        for l=1:Vxf0.L+1
            Vxf0.P(:,:,l) = eye(Vxf0.d);
        end
    end
    

    % Learn energy function
    t_start = tic;
    Vxf = learnEnergy(Vxf0,clfdm_data,options);

    % Learn GMM
    [Priors0, Mu0, Sigma0] = EM_init_kmeans(clfdm_data, n_gmm_components);
    [Priors_EM, Mu_EM, Sigma_EM, nbStep, loglik] = EM(clfdm_data, Priors0, Mu0, Sigma0);
    t_algo = toc(t_start);

    % Get unsafe DS handle
    in = 1:Vxf.d;
    out = Vxf.d+1:2*Vxf.d;
    fn_handle_GMR = @(x) GMR(Priors_EM, Mu_EM, Sigma_EM, x, in,out);
    % Get safe DS handle
    fn_handle = @(x) DS_stabilizer(x,fn_handle_GMR,Vxf,rho0,kappa0);

    
    x1_range = linspace(axis_bounds(1, 1), axis_bounds(1, 2), n_interpolation);
    x2_range = linspace(axis_bounds(2, 1), axis_bounds(2, 2), n_interpolation);

    % Create V contour
    [X, Y] = meshgrid(x1_range, x2_range);
    x = [X(:) Y(:)]';
    V = computeEnergy(x, [], Vxf);
    v_surf = reshape(V, n_interpolation, n_interpolation);

    % Create reproductions
    reproductions = {};
    initial_states = data_pos(:, data_start_idx);
    for i=1:n_demos
        curr_x = initial_states(:, i);
        curr_traj = [curr_x];

        for j=1:n_sim_steps
            dx = fn_handle(curr_x);
            curr_x = curr_x + dx*dt;
            curr_traj = [curr_traj curr_x];
        end
        reproductions{end+1} = curr_traj;
    end

    x1_dot_surf = zeros(n_interpolation);
    x2_dot_surf = zeros(n_interpolation);

    for i=1:n_interpolation
        for j=1:n_interpolation

            x_curr = [x1_range(i); x2_range(j)];
            x_dot_curr = fn_handle(x_curr);
            x1_dot_surf(j, i) = x_dot_curr(1);
            x2_dot_surf(j, i) = x_dot_curr(2);
        end
    end

    target_data = {};
    for i=1:n_demos
        if i==n_demos
            target_data{end+1} = data_pos(:, data_start_idx(i):end);
        else
            target_data{end+1} = data_pos(:, data_start_idx(i):data_start_idx(i+1));
        end

    end

    % Save parameters used for evaluation
    algoname = 'dmp';
    algo_params.L = Vxf0.L;
    algo_params.rho0 = rho0;
    algo_params.kappa0 = kappa0;
    algo_params.options = options;
    algo_params.n_gmm_components = n_gmm_components;
    algo_params.Vxf0 = Vxf0;
    algo_params.Vxf = Vxf;
    algo_params.Priors_EM = Priors_EM;
    algo_params.Mu_EM = Mu_EM;
    algo_params.Sigma_EM = Sigma_EM;
    algo_params.rand_init = rand_init;
    algo_params.seed = seed;


    data_params.n_samples = n_samples;
    data_params.n_demos = n_demos;

    output_file = fullfile(output_dir, strcat('plot_data_', algoname, '_', shapename, '.mat'));

    save(output_file, 'shapename', 'algoname', 'algo_params', 'data_params', ...
        'v_surf', 'x1_dot_surf', 'x2_dot_surf', 'reproductions', ...
        'x1_range', 'x2_range', 't_algo', 'target_data');

    fprintf("Seed %i: Finished shape %s \n", seed, shapename);

end






