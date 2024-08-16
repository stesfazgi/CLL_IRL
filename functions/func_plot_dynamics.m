function [x1_dot_surf, x2_dot_surf] = func_plot_dynamics(x, f_cl, axisBounds, n_interpolation, fig_title, savepath, data_pos)

    %%%%%%%%%%%% Default initializors
    if ~exist('n_interpolation', 'var')
        n_interpolation = 100;
    end
    
    
    x_range = linspace(min(axisBounds(:)), max(axisBounds(:)), n_interpolation);
    
    x1_dot_surf = zeros(n_interpolation);
    x2_dot_surf = zeros(n_interpolation);
    
    for i=1:n_interpolation
        for j=1:n_interpolation
            
            x_curr = [x_range(i); x_range(j)];
            
            x_dot_curr = replace(f_cl, x, x_curr);
            
            x1_dot_surf(j, i) = x_dot_curr(1);
            x2_dot_surf(j, i) = x_dot_curr(2);
            
        end
    end
    
    fig = figure('Position', [10 10 400 400]);
    hold on;
    if exist('data_pos', 'var')
        plot(data_pos(1, :), data_pos(2, :), 'black.');
    end
    l = streamslice(x_range, x_range, x1_dot_surf, x2_dot_surf);
    hold on;
    set(l, 'LineWidth', 1.2);
    set(l, 'Color', [.0, .0, .0]);
    xlabel('x1');
    ylabel('x2');
    xticks([-1 -0.5 0 0.5 1]);
    yticks([-1 -0.5 0 0.5 1]);
    xlim(axisBounds(1, :));
    ylim(axisBounds(2, :));
%     legend("Dataset", "Open Loop");
    if exist('fig_title', 'var')
        title(fig_title);
    end
    
    if exist('savepath', 'var')
        saveas(fig, savepath)
    end

end