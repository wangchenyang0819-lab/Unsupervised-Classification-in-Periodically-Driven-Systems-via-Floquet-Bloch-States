function combinedPlot()
    %% Part 1: Create left-side subplots (a, b, c, d)
    figure('Position', [100, 100, 1200, 800], 'Color', 'w');
    
    % Load data
    N = 13;
    try
        load('E:\mnist\kernelAIII11.mat', 'xx', 'yy', 'zz');
        u = xx; v = yy; w = zz;
        
        load('E:\mnist\kernelAIII31.mat', 'xx', 'yy', 'zz');
        uu = xx; vv = yy; ww = zz;
    catch
        % Generate sample data if files are not found
        [xx, yy, zz] = meshgrid(linspace(-1,1,N), linspace(-1,1,N), linspace(-1,1,1));
        u = xx; v = yy; w = sin(2*pi*xx).*cos(2*pi*yy);
        uu = sin(2*pi*xx); vv = cos(2*pi*yy); ww = xx.*yy;
    end
    
    [x, y, z] = meshgrid(0:2/(N-1):2, 0:2/(N-1):2, 1);
    
    % Subplot a: SSH model (dots + black arrows)
    ax_a = subplot('Position', [0.05, 0.55, 0.10, 0.40]);
    plotSSHVector(ax_a);
    text(ax_a, -0.45, 1.05, '(a)', 'Units', 'normalized', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    
    % Subplot b: Constant vector field (dots + black arrows)
    ax_b = subplot('Position', [0.05, 0.05, 0.10, 0.40]);
    plotConstantVector(ax_b);
    
    % Subplot c: Floquet model (colored dots + black arrows)
    ax_c = subplot('Position', [0.16, 0.55, 0.35, 0.40]);
    plotFloquetVector(ax_c, x, y, z, u, v, w);
    text(ax_c, -0.15, 1.05, '(b)', 'Units', 'normalized', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    
    % Subplot d: Another Floquet model (colored dots + black arrows)
    ax_d = subplot('Position', [0.16, 0.05, 0.35, 0.40]);
    plotFloquetVector(ax_d, x, y, z, uu, vv, ww);
    
    % Calculate positions for arrows and labels (currently unused)
    ax_c_pos = get(ax_c, 'Position');
    ax_d_pos = get(ax_d, 'Position');
    left_mid_x = ax_c_pos(1);
    left_mid_y = (ax_c_pos(2) + ax_d_pos(2) + ax_d_pos(4)) / 2;
    right_mid_x = ax_c_pos(1) + ax_c_pos(3);
    right_mid_y = (ax_c_pos(2) + ax_d_pos(2) + ax_d_pos(4)) / 2;
    mid_x = (left_mid_x + right_mid_x) / 2;
    mid_y = left_mid_y;
    
    %% Part 2: Compute data for right-side subplots
    [tpp, e, score, jiange] = computeSecondPart();
    
    %% Part 3: Create right-side subplots (e, f)
    % Subplot e: Phase surface
    ax_e = subplot('Position', [0.50, 0.55, 0.35, 0.40]);
    plotPhaseSurface(ax_e, tpp, jiange);
    
    % Subplot f: Eigenvalue plot with embedded PCA
    ax_f = subplot('Position', [0.50, 0.05, 0.35, 0.40]);
    plotEigenvaluesWithPCA(ax_f, e, score);
    
    %% Nested function definitions
    function plotSSHVector(ax)
        set(ax, 'NextPlot', 'add');
        y1 = linspace(0, 2, 11);
        x1 = zeros(size(y1));
        v11 = 0.8;
        w11 = 0.3;
        u1 = v11 + w11 * cos(pi * y1);
        v1 = w11 * sin(pi * y1);
        magnitudes = sqrt(u1.^2 + v1.^2);
        u1 = u1 ./ magnitudes;
        v1 = v1 ./ magnitudes;
        
        % Plot colored dots
        scatter(ax, x1, y1, 80, [0.25, 0.75, 0.95], 'filled', 'MarkerEdgeColor', 'k');
        
        % Plot black arrows
        quiver(ax, x1, y1, u1, v1, 0.25, 'LineWidth', 1.5, 'AutoScale', 'off', 'ShowArrowHead', 'on', 'Color', 'k');
        
        ylim(ax, [-0.1, 2.1]);
        axis(ax, 'equal');
        
        % Add gray horizontal reference line (k/π = 1)
        plot(ax, [min(x1), max(x1)], [1, 1], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
        
        set(ax, 'XTick', [], 'XTickLabel', []);
        set(ax, 'TickLabelInterpreter', 'latex', 'YTick', [0 1 2], 'FontSize', 28);
        grid(ax, 'on');
        set(ax, 'LineWidth', 3, 'Box', 'on');
        set(ax, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top');
        ylabel(ax, '$k/\pi$', 'Interpreter', 'latex', 'FontSize', 28, 'Color', 'k');
        ax.YAxis.TickLabelColor = 'k';
        ax.XAxis.TickLabelColor = 'k';
    end

    function plotConstantVector(ax)
        set(ax, 'NextPlot', 'add');
        y2 = linspace(0, 2, 11);
        x2 = zeros(size(y2));
        v11 = 0.3;
        w11 = 1.9;
        u2 = v11 + w11 * cos(pi * y2);
        v2 = w11 * sin(pi * y2);
        magnitudes = sqrt(u2.^2 + v2.^2);
        u2 = u2 ./ magnitudes;
        v2 = v2 ./ magnitudes;
        
        % Plot colored dots
        scatter(ax, x2, y2, 80, [0.25, 0.75, 0.95], 'filled', 'MarkerEdgeColor', 'k');
        
        % Plot black arrows
        quiver(ax, x2, y2, u2, v2, 0.25, 'LineWidth', 1.5, 'AutoScale', 'off', 'ShowArrowHead', 'on', 'Color', 'k');
        
        ylim(ax, [-0.1, 2.1]);
        axis(ax, 'equal');
        
        % Add gray horizontal reference line (k/π = 1)
        plot(ax, [min(x2), max(x2)], [1, 1], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
        
        set(ax, 'XTick', [], 'XTickLabel', []);
        set(ax, 'TickLabelInterpreter', 'latex', 'YTick', [0 1 2], 'FontSize', 28);
        grid(ax, 'on');
        set(ax, 'LineWidth', 3, 'Box', 'on');
        set(ax, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top');
        ylabel(ax, '$k/\pi$', 'Interpreter', 'latex', 'FontSize', 28, 'Color', 'k');
        ax.YAxis.TickLabelColor = 'k';
        ax.XAxis.TickLabelColor = 'k';
    end

    function plotFloquetVector(ax, x, y, z, u, v, w)
        set(ax, 'NextPlot', 'add');
        
        % Extract 2D positions
        X = x(:,:,1);
        Y = y(:,:,1);
        W = w(:,:,1);
        U = u(:,:,1);
        V = v(:,:,1);
        
        % Flatten vectors
        X_flat = X(:);
        Y_flat = Y(:);
        W_flat = W(:);
        U_flat = U(:);
        V_flat = V(:);
        
        % Create custom colormap
        cmap = createCustomColorMap();
        
        % Plot colored dots (z-component)
        scatter(ax, X_flat, Y_flat, 80, W_flat, 'filled', 'MarkerEdgeColor', 'k');
        colormap(ax, cmap);
        caxis(ax, [-1, 1]);
        
        % Plot black arrows
        quiver(ax, X_flat, Y_flat, U_flat/4, V_flat/4, 0.5, 'LineWidth', 1.5, 'AutoScale', 'off', 'ShowArrowHead', 'on', 'Color', 'k');
        
        % Add gray horizontal reference line (k/π = 1)
        plot(ax, [min(X_flat), max(X_flat)], [1, 1], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
        
        % Add colorbar
        c = colorbar(ax, 'eastoutside');
        c.Ticks = [-1, 1];
        c.TickLabels = {'$-1$', '$1$'};
        c.TickLabelInterpreter = 'latex';
        c.FontName = 'Times New Roman';
        c.FontSize = 26;
        c.Color = 'k';
        c.Label.Interpreter = 'latex';
        c.Label.FontSize = 26;
        c.Label.String = '$z$';
        c.Label.Position(1) = c.Label.Position(1) - 0.1;
        
        % Axis settings
        ylim(ax, [-0.1, 2.1]);
        axis(ax, 'equal');
        set(ax, 'XTick', [], 'XTickLabel', []);
        set(ax, 'TickLabelInterpreter', 'latex', 'YTick', [0 1 2], 'FontSize', 28);
        set(ax, 'LineWidth', 3, 'Box', 'on');
        set(ax, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top');
        grid(ax, 'on');
        ylabel(ax, '$k/\pi$', 'Interpreter', 'latex', 'FontSize', 28, 'Color', 'k');
        ax.YAxis.TickLabelColor = 'k';
        ax.XAxis.TickLabelColor = 'k';
        
        % Add time axis labels and arrows
        min_x = min(X_flat);
        max_x = max(X_flat);
        mid_x = (min_x + max_x) / 2;
        arrow1_x = (min_x + mid_x) / 2;
        arrow2_x = (mid_x + max_x) / 2;
        weiyi = -0.4;
        
        text(ax, min_x, weiyi, '$0$', 'Units', 'data', 'Interpreter', 'latex', 'FontSize', 26, 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
        text(ax, arrow1_x, weiyi, '$\longrightarrow$', 'Units', 'data', 'Interpreter', 'latex', 'FontSize', 26, 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
        text(ax, mid_x, weiyi, '$t$', 'Units', 'data', 'Interpreter', 'latex', 'FontSize', 26, 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
        text(ax, arrow2_x, weiyi, '$\longrightarrow$', 'Units', 'data', 'Interpreter', 'latex', 'FontSize', 26, 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
        text(ax, max_x, weiyi, '$T$', 'Units', 'data', 'Interpreter', 'latex', 'FontSize', 26, 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
    end
    
    function cmap = createCustomColorMap()
        % Create vertical gradient colormap
        positions = [0.0, 0.25, 0.5, 0.75, 1.0];
        colors = [0.00, 0.15, 0.65;  % Light mint green (z=1)
                  0.00, 0.50, 0.95;  % Light cyan
                  0.25, 0.75, 0.95;  % Blue-green
                  0.50, 0.95, 0.85;  % Turquoise
                  0.75, 1.00, 0.75]; % Deep navy blue (z=-1)
        n_colors = 256;
        cmap = interp1(positions, colors, linspace(0, 1, n_colors));
        cmap = min(1, max(0, cmap));
    end
    
    function [tpp, e, score, jiange] = computeSecondPart()
        load('E:\mnist\AIII.mat','e','score','tpp','B','f');
        jiange = 600;  % Resolution for phase diagram
    end
    
    function plotPhaseSurface(ax, tpp, jiange)
        tick_positions = [1, jiange/2, 0.95*jiange];
        color1 = [194/256, 206/256, 220/256];  % Light red (region 1)
        color2 = [145/256, 173/256, 158/256];  % Light green (region 2)
        color3 = [216/256, 156/256, 122/256];  % Blue (region 3)
        color4 = [0.5, 0, 0.5];                % Dark purple (region 4)
        
        % Create RGB matrix
        R = zeros(size(tpp));
        G = zeros(size(tpp));
        B = zeros(size(tpp));
        R(tpp == 1) = color1(1); G(tpp == 1) = color1(2); B(tpp == 1) = color1(3);
        R(tpp == 2) = color2(1); G(tpp == 2) = color2(2); B(tpp == 2) = color2(3);
        R(tpp == 3) = color3(1); G(tpp == 3) = color3(2); B(tpp == 3) = color3(3);
        R(tpp == 4) = color4(1); G(tpp == 4) = color4(2); B(tpp == 4) = color4(3);
        
        RGB = cat(3, R, G, B);
        load('E:\mnist\AIII.mat','B','f');
        hold on;
        image(ax, RGB);
        hold on;
        
        % Axis settings
        set(ax, 'FontName', 'Times New Roman', 'FontSize', 18);
        axis(ax, 'square');
        set(ax, 'LineWidth', 2, 'Box', 'on');
        set(ax, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top');
        set(ax, 'TickLabelInterpreter', 'latex', 'XTick', tick_positions, 'XTickLabel', {'0', '0.5', '1'}, ...
             'YTick', tick_positions, 'YTickLabel', {'-0.5','0', '0.5'}, 'FontSize', 28);
        text(ax, -0.15, 1.05, '(c)', 'Units', 'normalized', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
        xlabel(ax, '$\theta_{\mathrm{Re}}T/\pi$', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k');
        ylabel(ax, '$\theta_{\mathrm{Im}} T/\pi$', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k');
        ax.YAxis.TickLabelColor = 'k';
        ax.XAxis.TickLabelColor = 'k';
    end
    
    function plotEigenvaluesWithPCA(ax, e, score)
        set(ax, 'NextPlot', 'add');
        text(ax, -0.15, 1.05, '(d)', 'Units', 'normalized', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
        
        % Plot eigenvalues
        axis(ax, 'square');
        plot(ax, 1:8, e(1:8), 'o', 'Color', [0.25, 0.75, 0.95], 'MarkerFaceColor', [0.25, 0.75, 0.95], 'LineWidth', 1.5, 'MarkerSize', 8);
        ylim(ax, [0, 1.1]);
        set(ax, 'XTick', 1:8, 'XTickLabel', {'1','2','3','4','5','6','7','8'}, 'FontSize', 18, 'FontName', 'Times New Roman');
        set(ax, 'TickLabelInterpreter', 'latex', 'YTick', [0 1], 'FontSize', 28);
        xlabel(ax, '$n$', 'Interpreter', 'latex', 'FontSize', 20);
        ylabel(ax, '$\lambda_n$', 'Interpreter', 'latex', 'FontSize', 20);
        box(ax, 'on');
        
        % Embedded PCA plot
        ax_pca = axes('Position', [0.58, 0.10, 0.12, 0.12]);
        colors = zeros(size(score, 1), 3);
        
        % Assign colors based on PCA scores
        for i = 1:size(score, 1)
            if score(i,1) <= -0.01
                colors(i, :) = [194/256, 206/256, 220/256];
            elseif score(i,1) >= -0.01 && score(i,2) <= 0
                colors(i, :) = [216/256, 156/256, 122/256];
            elseif score(i,1) >= -0.01 && score(i,2) >= -0
                colors(i, :) = [145/256, 173/256, 158/256];
            elseif score(i,1) > 0.07
                colors(i, :) = [0.5, 0, 0.5];
            else
                colors(i, :) = [0, 0, 0];
            end
        end
        
        hold(ax_pca, 'on');
        for i = 1:size(score, 1)
            scatter(ax_pca, score(i,1), score(i,2), 40, colors(i, :), 'filled');
        end
        hold(ax_pca, 'off');
        
        % PCA plot settings
        xlim(ax_pca, [-0.05, 0.1]);
        ylim(ax_pca, [-0.05, 0.1]);
        xlabel(ax_pca, 'PC1', 'FontSize', 10);
        ylabel(ax_pca, 'PC2', 'FontSize', 10);
        title(ax_pca, 'PCA Embedding', 'FontSize', 10);
        box(ax_pca, 'on');
        axis(ax_pca, 'equal');
    end
end