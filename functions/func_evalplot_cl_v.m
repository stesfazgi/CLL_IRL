function [] = func_evalplot_cl_v(x1_range, x2_range, v_surf, x1dot, x2dot, reproductions, target_data, algoname, filename)

    pdf_plt_size = 7.5;

    
    %%%%%%%%%%%%%%%%%% PLOT STREAMLINE AND REPRODUCTIONS %%%%%%%%%%%%%%%%%%
    
    fig = figure('Position', [0 0 400 400], 'Visible',"off");
    hold on;
    l = streamslice(x1_range, x2_range, x1dot, x2dot);

    hold on;
    
    for i=1:length(target_data)
        curr_traj = target_data{i};
        plot(curr_traj(1, :), curr_traj(2, :), "black.", 'LineWidth', 2);
        hold on;
    end
    
    if ~strcmp(algoname, 'ol')
        for i=1:length(reproductions)
            curr_traj = reproductions{i};
            plot(curr_traj(1, :), curr_traj(2, :), "r", 'LineWidth', 2);
            hold on;
        end
    end
    

    set(l, 'LineWidth', 1.2);
    set(l, 'Color', [.4, .4, .4]);
    xlim([min(x1_range), max(x1_range)]);
    ylim([min(x2_range), max(x2_range)]);
    xticks([]);
    yticks([]);
    
    box on;
    

    set(gcf, 'PaperPosition', [0 0 pdf_plt_size pdf_plt_size]);
    set(gcf, 'PaperSize', [pdf_plt_size pdf_plt_size]);
    set(gca,'position',[0 0 1 1],'units','normalized');
    set(gca, 'LineWidth', 2);
    cl_plot_filename = strcat(filename, "_CL", ".pdf");
    saveas(gcf, cl_plot_filename, 'pdf')
    
    if strcmp(algoname, 'ol')
        return;
    end
    
    %%%%%%%%%%%%%%%%%% PLOT CONTOUR OF V %%%%%%%%%%%%%%%%%%
    
    % Find largest value of V for all trajectory starts
    largest_V = 0;
    smallest_V = 1e10;
    for i=1:length(reproductions)
        curr_traj = squeeze(reproductions{i});
        [~, idx1] = min(abs(x1_range - curr_traj(1, 1)));
        [~, idx2] = min(abs(x2_range - curr_traj(2, 1)));

        if v_surf(idx2, idx1) > largest_V
            largest_V = v_surf(idx2, idx1);
        end

        [~, idx1] = min(abs(x1_range - curr_traj(1, end)));
        [~, idx2] = min(abs(x2_range - curr_traj(2, end)));

        if v_surf(idx2, idx1) < smallest_V
            smallest_V = v_surf(idx2, idx1);
        end
    end
    
    v_surf_part1 = v_surf .* (v_surf < largest_V);
    v_surf_part2 = v_surf .* (v_surf >= largest_V);

    % Shift both parts so that the smallest value is 0
    v_surf_part1 = v_surf_part1 - min(v_surf_part1(:));
    % v_surf_part2 = v_surf_part2 - largest_V;
    % Normalize both parts
    v_surf_part1 = v_surf_part1 ./ max(v_surf_part1(:));
    % v_surf_part2 = v_surf_part2 ./ max(v_surf_part2(:));
    % Apply mask to part1
    % v_surf_part1 = v_surf_part1 + 1*(v_log >= largest_V);

    v_surf_final = v_surf_part1 + 1*(v_surf_part2 > 0);

    % PART 1
    n_inter_B_G = 9;
    n_inter_G_Y = 3;
    % PART 2
    n_inter_Y_R = 9;

    interp_cmap_B_G = [zeros(n_inter_B_G, 1), linspace(0, 1, n_inter_B_G)', linspace(1, 0, n_inter_B_G)'];
    interp_cmap_G_Y = [linspace(0, 1, n_inter_G_Y)', ones(n_inter_G_Y, 1), zeros(n_inter_G_Y, 1)];

    interp_cmap_Y_R = [ones(n_inter_Y_R, 1), linspace(1, 0, n_inter_Y_R)', zeros(n_inter_Y_R, 1)];

    cmap_final = ...
        [interp_cmap_B_G;
         interp_cmap_G_Y;
         interp_cmap_Y_R];

    fig = figure('Position', [0 0 400 400], 'Visible',"off");
    hold on;

    if strcmp(algoname, 'irl')
        s = pcolor(v_surf);
    else
        s = pcolor(v_surf_final);
    end
    
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    colormap(cmap_final);

%     if strcmp(algoname, 'irl')
%         contour_levels = [linspace(min(v_surf(:)), max(v_surf(:)), 100)];
%     else
%         contour_levels = linspace(smallest_V, largest_V, 50);
%         contour_levels = [contour_levels logspace(log10(largest_V), log10(max(v_surf(:))), 20)];
%     end
% 
% 
%     contour(x1_range, x2_range, v_surf, contour_levels, "b");


    

    hold on;
    box on;
    
    xticks([]);
    yticks([]);
    
    set(gcf, 'PaperPosition', [0 0 pdf_plt_size pdf_plt_size]);
    set(gcf, 'PaperSize', [pdf_plt_size pdf_plt_size]);
    set(gca,'position',[0 0 1 1],'units','normalized');
    set(gca, 'LineWidth', 2);
    V_plot_filename = strcat(filename, "_V", ".pdf");
    saveas(gcf, V_plot_filename, 'pdf')

end

