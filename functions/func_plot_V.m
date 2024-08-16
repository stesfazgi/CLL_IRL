function [outputArg1,outputArg2] = func_plot_V(x, V_star, axisBounds, figtitle, data_pos)

    if ~exist('data_pos', 'var')
        data_pos = [];
    end

    n_interpolation = 100;
    n_contour_levels = 100;
    plot_width = 600;

    x_range = linspace(min(axisBounds(:)), max(axisBounds(:)), n_interpolation);
    v_surf = zeros(n_interpolation);

    % Evaluate V on n_interpolation by n_interpolation grid
    for i=1:n_interpolation
        for j=1:n_interpolation
            
            x_curr = [x_range(i); x_range(j)];
            
            v_surf(j, i) = replace(V_star, x, x_curr);
            
        end
    end

    
    fig = figure('Position', [10 10 plot_width 600], 'Visible', 'on');
    clfsos_levels = logspace(-2, log10(max(v_surf(:))), n_contour_levels);
    contour(x_range, x_range, v_surf, clfsos_levels, 'HandleVisibility', 'off');
    hold on;
    xlabel('x1');
    ylabel('x2');
    xlim(axisBounds(1, :));
    ylim(axisBounds(2, :));
    
    if size(data_pos, 2) > 0
        plot(data_pos(1, :), data_pos(2, :), 'r+');
    end
    
    if exist('figtitle', 'var')
        title(figtitle);
    end
    
end

