function combinedPlot()
    %% Create figure with 4 subplots (2x2 grid)
    figure('Position', [100, 100, 950, 850], 'Color', 'w', 'Name', 'Figure 17', 'NumberTitle', 'off');
    
    % Load data for eigenvalue plot
    load('E:\mnist\A1.mat', 'e', 'score');
    
    % Subplot (a): Lattice model (top-left)
    ax_a = subplot(2, 2, 1);
    plotLatticeModel(ax_a, e, score);
    text(ax_a, -0.15, 1.05, '(a)', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    
    % Subplot (b): Phase diagram (top-right)
    ax_b = subplot(2, 2, 2);
    plotPhaseDiagram(ax_b, score, e);
    text(ax_b, -0.15, 1.05, '(b)', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    
    % Subplot (c): Trivial link (bottom-left)
    ax_c = subplot(2, 2, 3);
    plotTrivialLink(ax_c);
    text(ax_c, -0.15, 1.05, '(c)', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    
    % Subplot (d): Non-trivial link (bottom-right)
    ax_d = subplot(2, 2, 4);
    plotNonTrivialLink(ax_d);
    text(ax_d, -0.15, 1.05, '(d)', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    
    %% ===== Nested Function Definitions =====
    
    function plotLatticeModel(ax, e, score)
        % Plot lattice model with couplings
        
        hold(ax, 'on');
        axis(ax, 'equal');
        set(ax, 'XTick', [], 'YTick', [], 'Box', 'on');
        
        % Create 4x4 coordinate grid
        [x_grid, y_grid] = meshgrid(-2:1, -2:1);
        
        % Draw grid lines
        for i = -2:1
            plot(ax, [-2, 1], [i, i], 'k-', 'LineWidth', 1.5);  % Horizontal lines
            plot(ax, [i, i], [-2, 1], 'k-', 'LineWidth', 1.5);  % Vertical lines
        end
        
        % Draw lattice points
        even_mask = mod(x_grid + y_grid, 2) == 0;
        odd_mask = mod(x_grid + y_grid, 2) == 1;
        
        scatter(ax, x_grid(even_mask), y_grid(even_mask), 120, 'k', 'filled', 'MarkerEdgeColor', 'k');
        scatter(ax, x_grid(odd_mask), y_grid(odd_mask), 120, 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        
        % Define colors for different couplings
        colors = struct('t1', [1 0 0], ...    % Red
                        't2', [0 0.7 0], ...  % Green
                        't3', [0 0 1], ...    % Blue
                        't4', [1 0.5 0]);     % Orange
        
        % Draw couplings and labels
        for i = 1:numel(x_grid)
            if odd_mask(i)
                xi = x_grid(i);
                yi = y_grid(i);
                
                % Right coupling (t1, red)
                if xi < 1 && any(x_grid(even_mask) == xi+1 & y_grid(even_mask) == yi)
                    text(ax, xi + 0.5, yi - 0.15, '$J$', 'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle', 'FontSize', 16, 'Color', 'k', ...
                         'Interpreter', 'latex', 'FontWeight', 'bold');
                    plot(ax, [xi, xi+1], [yi, yi], 'Color', colors.t1, 'LineWidth', 4);
                end
                
                % Down coupling (t4, orange)
                if yi > -2 && any(x_grid(even_mask) == xi & y_grid(even_mask) == yi-1)
                    text(ax, xi + 0.15, yi - 0.5, '$J$', 'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle', 'FontSize', 16, 'Color', 'k', ...
                         'Interpreter', 'latex', 'FontWeight', 'bold');
                    plot(ax, [xi, xi], [yi, yi-1], 'Color', colors.t4, 'LineWidth', 4);
                end
                
                % Left coupling (t3, blue)
                if xi > -2 && any(x_grid(even_mask) == xi-1 & y_grid(even_mask) == yi)
                    text(ax, xi - 0.5, yi + 0.15, '$J$', 'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle', 'FontSize', 16, 'Color', 'k', ...
                         'Interpreter', 'latex', 'FontWeight', 'bold');
                    plot(ax, [xi, xi-1], [yi, yi], 'Color', colors.t3, 'LineWidth', 4);
                end
                
                % Up coupling (t2, green)
                if yi < 1 && any(x_grid(even_mask) == xi & y_grid(even_mask) == yi+1)
                    text(ax, xi - 0.15, yi + 0.5, '$J$', 'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle', 'FontSize', 16, 'Color', 'k', ...
                         'Interpreter', 'latex', 'FontWeight', 'bold');
                    plot(ax, [xi, xi], [yi, yi+1], 'Color', colors.t2, 'LineWidth', 4);
                end
            end
        end
        
        % Set axis limits and hide axes
        xlim(ax, [-2.2, 1.2]);
        ylim(ax, [-2.2, 1.2]);
        set(ax, 'Box', 'off');
        set(ax, 'XColor', 'none', 'YColor', 'none');
    end
    
    function plotPhaseDiagram(ax, score, e)
        % Plot phase diagram with three regions
        
        set(ax, 'Box', 'on');
        hold(ax, 'on');
        axis(ax, 'square');
        
        % Define axis ranges
        x_min = 0;
        x_max = 3;
        y_min = -1;
        y_max = 1.0;
        
        % Define region boundaries
        x1 = 1.1;
        x2 = 2.0;
        
        % Define colors for different regions
        color1 = [194/256, 206/256, 220/256];  % Light blue (trivial phase)
        color2 = [145/256, 173/256, 158/256];  % Light green (π phase)
        color3 = [216/256, 156/256, 122/256];  % Orange (0π phase)
        color4 = [0.2, 0.2, 0.4];              % Dark blue for PC1
        color5 = [0.4, 0.2, 0.2];              % Dark red for PC2
        white_color = [1, 1, 1];               % White
        
        % Plot three regions with gradient
        patch(ax, 'XData', [x_min, x_min, x1, x1], ...
              'YData', [y_min, y_max, y_max, y_min], ...
              'FaceVertexCData', [color1; white_color; white_color; color1], ...
              'FaceColor', 'interp', 'EdgeColor', 'none');
        
        patch(ax, 'XData', [x1, x1, x2, x2], ...
              'YData', [y_min, y_max, y_max, y_min], ...
              'FaceVertexCData', [color2; white_color; white_color; color2], ...
              'FaceColor', 'interp', 'EdgeColor', 'none');
        
        patch(ax, 'XData', [x2, x2, x_max, x_max], ...
              'YData', [y_min, y_max, y_max, y_min], ...
              'FaceVertexCData', [color3; white_color; white_color; color3], ...
              'FaceColor', 'interp', 'EdgeColor', 'none');
        
        % Add labels for each phase
        x_center1 = (x_min + x1) / 2;
        x_center2 = (x1 + x2) / 2;
        x_center3 = (x2 + x_max) / 2;
        vertical_position = 0.7;
        font_size = 14;
        text_color = [0, 0, 0];
        
        text(ax, x_center1, vertical_position, 'trivial phase', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', font_size, 'Color', text_color, ...
             'Interpreter', 'latex', 'FontWeight', 'bold');
        
        text(ax, x_center2, vertical_position, '$\pi$ phase', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', font_size, 'Color', text_color, ...
             'Interpreter', 'latex', 'FontWeight', 'bold');
        
        text(ax, x_center3, vertical_position, '$0\pi$ phase', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', font_size, 'Color', text_color, ...
             'Interpreter', 'latex', 'FontWeight', 'bold');
        
        % Plot PCA components
        pc1 = score(:, 1);
        pc2 = score(:, 2);
        x_vals = 3/31:3/31:90/31;
        
        plot(ax, x_vals, 0.3*pc1./max(pc1) + 0.1, 'o', 'Color', color4, ...
             'MarkerSize', 6, 'LineStyle', 'none', 'LineWidth', 1.5, 'MarkerFaceColor', color4);
        hold on
        plot(ax, x_vals, 0.3*pc2./max(pc2) + 0.2, 's', 'Color', color5, ...
             'MarkerSize', 6, 'LineStyle', 'none', 'LineWidth', 1.5, 'MarkerFaceColor', color5);
        
        % Set axis properties
        set(ax, 'LineWidth', 2, 'Box', 'on');
        set(ax, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top');
        ax.XAxis.TickLabelColor = 'k';
        ax.YAxis.TickLabelColor = 'k';
        
        xlim(ax, [x_min, x_max]);
        ylim(ax, [y_min, y_max]);
        set(ax, 'XTick', [0, x1, x2, x_max], 'FontSize', 18);
        set(ax, 'TickLabelInterpreter', 'latex', 'XTickLabel', {'', num2str(x1), num2str(x2), num2str(x_max)}, 'FontSize', 18);
        set(ax, 'YTick', [], 'FontSize', 18);
        xlabel(ax, '$JT/\pi$', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'k');
        grid(ax, 'on');
        hold(ax, 'off');
        
        % Create inset for eigenvalues
        ax_e = axes('Position', [0.63, 0.63, 0.10, 0.10]);
        plot(ax_e, 1:8, abs(e(1:8)), 'o', 'Color', [0.25, 0.75, 0.95], ...
            'MarkerFaceColor', [0.25, 0.75, 0.95], ...
            'LineWidth', 0.5, 'MarkerSize', 4);
        ylim(ax_e, [0, 1.1]);
        xlabel(ax_e, '$n$', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 10);
        ylabel(ax_e, '$\lambda_n$', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 10);
        box(ax_e, 'on');
    end
    
    function plotTrivialLink(ax)
        % Plot trivial link diagram
        
        % Load trivial link data
        load('E:\mnist\trivallink.mat', 'redPoints', 'N', 'bluePoints');
        
        set(ax, 'NextPlot', 'add');
        set(ax, 'FontName', 'Times New Roman', 'FontSize', 16);
        
        % Plot red points
        if ~isempty(redPoints)
            plot(ax, redPoints(:,1)/N*2-1, redPoints(:,2)/N*2-1, 'r.', 'MarkerSize', 8);
        end
        
        % Plot blue points
        if ~isempty(bluePoints)
            plot(ax, bluePoints(:,1)/N*2-1, bluePoints(:,2)/N*2-1, 'b.', 'MarkerSize', 8);
        end
        
        % Calculate vertices of the quadrilateral
        A = [1, 1; 1, -1];
        b1 = [1; 1];
        b2 = [1; -1];
        b3 = [-1; -1];
        b4 = [-1; 1];
        vertex1 = A\b1;  % (1,0)
        vertex2 = A\b2;  % (0,1)
        vertex3 = A\b3;  % (-1,0)
        vertex4 = A\b4;  % (0,-1)
        
        % Draw dashed quadrilateral
        vertices = [vertex1, vertex2, vertex3, vertex4, vertex1];
        plot(ax, vertices(1,:), vertices(2,:), 'k--', 'LineWidth', 1.5);
        
        % Set axis properties
        axis(ax, 'square');
        set(ax, 'LineWidth', 2, 'Box', 'on');
        set(ax, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top');
        ax.XAxis.TickLabelColor = 'k';
        ax.YAxis.TickLabelColor = 'k';
        xlabel(ax, '$k_x/\pi$', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 20);
        ylabel(ax, '$k_y/\pi$', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 20);
        set(ax, 'TickLabelInterpreter', 'latex', 'XTick', [-0.90 0 0.90], 'XTicklabel', [-1 0 1], 'FontSize', 18);
        set(ax, 'TickLabelInterpreter', 'latex', 'YTick', [-0.90 0 0.90], 'YTicklabel', [-1 0 1], 'FontSize', 18);
    end
    
    function plotNonTrivialLink(ax)
        % Plot non-trivial link diagram
        
        % Load non-trivial link data
        load('E:\mnist\nontrivallink.mat', 'redPoints', 'N', 'bluePoints');
        
        set(ax, 'NextPlot', 'add');
        set(ax, 'FontName', 'Times New Roman', 'FontSize', 16);
        
        % Plot blue points in certain regions
        for i = 1:size(bluePoints, 1)
            plot(ax, bluePoints(i,1)/N*2-1, bluePoints(i,2)/N*2-1, 'b.', 'MarkerSize', 8);
        end
        
        % Plot red points
        for i = 1:size(redPoints, 1)
            plot(ax, redPoints(i,1)/N*2-1, redPoints(i,2)/N*2-1, 'r.', 'MarkerSize', 8);
        end
        
        % Additional blue points in specific x ranges
        for i = 1:size(bluePoints, 1)
            if bluePoints(i,1)/N*2-1 > -0.29 && bluePoints(i,1)/N*2-1 < -0.2
                plot(ax, bluePoints(i,1)/N*2-1, bluePoints(i,2)/N*2-1, 'b.', 'MarkerSize', 8);
            end
        end
        for i = 1:size(bluePoints, 1)
            if bluePoints(i,1)/N*2-1 > 0.7 && bluePoints(i,1)/N*2-1 < 0.82
                plot(ax, bluePoints(i,1)/N*2-1, bluePoints(i,2)/N*2-1, 'b.', 'MarkerSize', 8);
            end
        end
        
        % Calculate vertices of the quadrilateral
        A = [1, 1; 1, -1];
        b1 = [1; 1];
        b2 = [1; -1];
        b3 = [-1; -1];
        b4 = [-1; 1];
        vertex1 = A\b1;  % (1,0)
        vertex2 = A\b2;  % (0,1)
        vertex3 = A\b3;  % (-1,0)
        vertex4 = A\b4;  % (0,-1)
        
        % Draw dashed quadrilateral
        vertices = [vertex1, vertex2, vertex3, vertex4, vertex1];
        plot(ax, vertices(1,:), vertices(2,:), 'k--', 'LineWidth', 1.5);
        
        % Set axis properties
        axis(ax, 'square');
        set(ax, 'LineWidth', 2, 'Box', 'on');
        set(ax, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top');
        ax.XAxis.TickLabelColor = 'k';
        ax.YAxis.TickLabelColor = 'k';
        xlabel(ax, '$k_x/\pi$', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 20);
        ylabel(ax, '$k_y/\pi$', 'Interpreter', 'latex', 'Color', 'k', 'FontSize', 20);
        set(ax, 'TickLabelInterpreter', 'latex', 'XTick', [-0.90 0 0.90], 'XTicklabel', [-1 0 1], 'FontSize', 18);
        set(ax, 'TickLabelInterpreter', 'latex', 'YTick', [-0.90 0 0.90], 'YTicklabel', [-1 0 1], 'FontSize', 18);
    end
end