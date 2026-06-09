function combinedPlot()
    %% Part 1: Create left-side subplots (a, b)
    figure('Position', [100, 100, 1200, 800], 'Color', 'w');
    
    % Load data
    N = 13;
    
    load('E:\mnist\kernelD0fu1.mat', 'xx', 'yy', 'zz');
    u = xx; v = yy; w = zz;
    
    load('E:\mnist\kernelD151.mat', 'xx', 'yy', 'zz');
    uu = xx; vv = yy; ww = zz;
    
    % Adjust specific points
    uu(5,10) = -u(5,10);
    vv(5,10) = -v(5,10);
    ww(5,10) = -w(5,10);
    uu(14-5,10) = uu(5,10);
    vv(14-5,10) = -vv(5,10);
    ww(14-5,10) = -ww(5,10);
    
    [x, y, z] = meshgrid(0:2/(N-1):2, 0:2/(N-1):2, 1);
    
    % Subplot (a): Floquet model
    ax_a = subplot('Position', [0.05, 0.55, 0.35, 0.35]);
    plotFloquetVector(ax_a, x, y, z, u, v, w, false); % 不显示t/T标签
    text(ax_a, -0.2, 1.05, '(a)', 'Units', 'normalized', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    
    % Subplot (c): Another Floquet model
    ax_b = subplot('Position', [0.05, 0.1, 0.35, 0.35]);
    plotFloquetVector(ax_b, x, y, z, uu, vv, ww, false); % 不显示t/T标签
    text(ax_b, -0.2, 1.05, '(c)', 'Units', 'normalized', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    
    %% Part 2: Load data for right-side subplots
    load('E:\mnist\D.mat', 'e', 'score', 'tpp', 'TT');
    jiange = 600;
    
    %% Part 3: Create right-side subplots (b, d)
    % Subplot (b): Phase surface
    ax_c = subplot('Position', [0.40, 0.55, 0.35, 0.35]);
    plotPhaseSurface(ax_c, tpp, jiange);
    
    % Subplot (d): Eigenvalue diagram with PCA embedding
    ax_d = subplot('Position', [0.40, 0.1, 0.35, 0.35]);
    plotEigenvaluesWithPCA(ax_d, e, score);
    
    %% ===== Nested Function Definitions =====
    function plotFloquetVector(ax, x, y, z, u, v, w, showXLabel)
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
        cmap = [0.00, 0.15, 0.65;   % Light mint green
                0.00, 0.50, 0.95;  % Light cyan
                0.25, 0.75, 0.95;  % Blue-green
                0.50, 0.95, 0.85;  % Turquoise
                0.75, 1.00, 0.75]; % Deep navy blue
        
        % Plot colored dots representing z-component
        scatter(ax, X_flat, Y_flat, 80, W_flat, 'filled', 'MarkerEdgeColor', 'k');
        colormap(ax, cmap);
        caxis(ax, [-1, 1]);
        
        % Add black arrows
        quiver(ax, X_flat, Y_flat, U_flat/4, V_flat/4, 0.5, ...
            'LineWidth', 1.5, 'AutoScale', 'off', 'Color', 'k');
        
        % Add gray horizontal reference line
        plot(ax, [min(X_flat), max(X_flat)], [1, 1], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
        
        % Add colorbar
        c = colorbar(ax, 'eastoutside');
        c.Ticks = [-1, 1];
        c.TickLabels = {'$-1$', '$1$'};
        c.TickLabelInterpreter = 'latex';
        c.FontName = 'Times New Roman';
        c.FontSize = 26;
        c.Color = 'k';

        
        % Axis settings
        set(ax, 'XTick', [0 1 2], 'YTick', [0 1 2], 'FontSize', 28);
        set(ax, 'TickLabelInterpreter', 'latex', 'LineWidth', 3, 'Box', 'on');
        set(ax, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top');
        ylim(ax, [-0.1, 2.1]);
        axis(ax, 'square');
        grid(ax, 'on');
        
        ax.YAxis.TickLabelColor = 'k';
        ax.XAxis.TickLabelColor = 'k';
        
        % Y-axis label
        ylabel(ax, '$k/\pi$', 'Interpreter', 'latex', 'FontSize', 28, 'Color', 'k');
        
        % 删除t/T标签，只保留刻度标签
        set(ax, 'XTickLabel', {'$0$', '$T/2$', '$T$'});
        xticks = get(ax, 'XTick');
xticklabels = get(ax, 'XTickLabel');
for i = 1:length(xticks)
    text(ax, xticks(i), -0.3, xticklabels{i}, ...  % 调整y坐标
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
         'FontSize', 28, 'Interpreter', 'latex');
end
% 清除原来的刻度标签
set(ax, 'XTickLabel', {});
        % 删除xlabel调用
        % 完全删除xlabel代码行
    end
    
    function plotPhaseSurface(ax, tpp, jiange)
        % Define region colors
        colors = {[194/256, 206/256, 220/256];  % Region 1
                  [145/256, 173/256, 158/256];  % Region 2
                  [216/256, 156/256, 122/256];  % Region 3
                  [0.5, 0, 0.5]};               % Region 4
        
        % Create RGB matrix
        R = zeros(size(tpp));
        G = zeros(size(tpp));
        B = zeros(size(tpp));
        
        for i = 1:4
            idx = tpp == i;
            R(idx) = colors{i}(1);
            G(idx) = colors{i}(2);
            B(idx) = colors{i}(3);
        end
        
        % Display image
        RGB = cat(3, R, G, B);
        image(ax, RGB);
        hold on;
        
        % Axis settings
        tick_positions = [2, 2000, 3900];
        set(ax, 'FontName', 'Times New Roman', 'FontSize', 18, 'LineWidth', 2, 'Box', 'on');
        set(ax, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top');
        set(ax, 'TickLabelInterpreter', 'latex', ...
            'XTick', tick_positions, 'XTickLabel', {'$-2$', '$0$', '$2$'}, ...
            'YTick', tick_positions, 'YTickLabel', {'$2$', '$0$', '$-2$'}, ...
            'FontSize', 24);

        text(ax, -0.22, 1.05, '(b)', 'Units', 'normalized', 'FontSize', 20, ...
            'FontWeight', 'bold', 'FontName', 'Times New Roman');
        
        % Axis labels for phase surface
        xl=xlabel(ax, '$J_{1}T/\pi$', 'Interpreter', 'latex', 'FontSize', 28, 'Color', 'k');
        xl.Position(2) = xl.Position(2) -0.03*xl.Position(2);

        ylabel(ax, '$J_{2}T/\pi$', 'Interpreter', 'latex', 'FontSize', 28, 'Color', 'k');
        axis(ax, 'square');
        
        ax.YAxis.TickLabelColor = 'k';
        ax.XAxis.TickLabelColor = 'k';
    end
    
    function plotEigenvaluesWithPCA(ax, e, score)
        set(ax, 'NextPlot', 'add');
        text(ax, -0.22, 1.05, '(d)', 'Units', 'normalized', 'FontSize', 20, ...
            'FontWeight', 'bold', 'FontName', 'Times New Roman');
        
        % Plot eigenvalues
        axis(ax, 'square');
        plot(ax, 1:12, e(1:12), 'o', 'Color', [0.25, 0.75, 0.95], ...
            'MarkerFaceColor', [0.25, 0.75, 0.95], ...
            'LineWidth', 1.5, 'MarkerSize', 8);
        
        ylim(ax, [0, 1.1]);
        set(ax, 'XTick', 1:7:8, 'XTickLabel', {'$1$', '$8$'}, 'FontSize', 24, 'FontName', 'Times New Roman');
        set(ax, 'TickLabelInterpreter', 'latex', 'YTick', [0 0.2 0.4 0.6 0.8 1], 'FontSize', 24);
        
        % Axis labels for eigenvalue plot
        xlabel(ax, '$n$', 'Interpreter', 'latex', 'FontSize', 28, 'Color', 'k');
xlab = get(ax, 'XLabel');
current_pos = get(xlab, 'Position');
set(xlab, 'Position', [current_pos(1), current_pos(2) + 0.05, current_pos(3)]);
        ylabel(ax, '$\lambda_n$', 'Interpreter', 'latex', 'FontSize', 28, 'Color', 'k');
        box(ax, 'on');
        
        % Embedded PCA subplot
        ax_pca = axes('Position', [0.53, 0.20, 0.12, 0.12]);
        colors = zeros(size(score, 1), 3);
        
        % Define color mapping
        color_map = containers.Map({1,2,3,4,5,6,7,8}, ...
            {[194/256, 206/256, 220/256], ...  % Region 1
             [145/256, 173/256, 158/256], ...  % Region 2
             [145/256, 173/256, 158/256], ...  % Region 2
             [0.5, 0, 0.5], ...               % Region 4
             [216/256, 156/256, 122/256], ...  % Region 3
             [216/256, 156/256, 122/256], ...  % Region 3
             [194/256, 206/256, 220/256], ...  % Region 1
             [0.5, 0, 0.5]});                  % Region 4
        
        % Assign colors and plot PCA
        hold(ax_pca, 'on');
        for i = 1:size(score, 1)
            colors(i, :) = color_map(TT(i,1));
            
            if ismember(TT(i,1), [1,2,4,5])
                marker = '^';
            else
                marker = 'o';
            end
            scatter(ax_pca, score(i,1), score(i,2), 40, colors(i, :), 'filled', 'Marker', marker);
        end
        hold(ax_pca, 'off');
        
        % PCA plot settings
        xlim(ax_pca, [-0.1, 0.15]);
        ylim(ax_pca, [-0.1, 0.15]);
        
        % Axis labels for PCA subplot
        xlabel(ax_pca, 'PC1', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(ax_pca, 'PC2', 'FontSize', 12, 'FontWeight', 'bold');
        title(ax_pca, 'PCA Embedding', 'FontSize', 12, 'FontWeight', 'bold');
        box(ax_pca, 'on');
        axis(ax_pca, 'equal');
        
        % Tick labels for PCA subplot
        set(ax_pca, 'XTick', [-0.1, 0, 0.1], 'XTickLabel', {'-0.1', '0', '0.1'}, 'FontSize', 10);
        set(ax_pca, 'YTick', [-0.1, 0, 0.1], 'YTickLabel', {'-0.1', '0', '0.1'}, 'FontSize', 10);
    end
end