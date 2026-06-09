clear

% Part 1: Create phase diagram
figure1 = figure;
axes1 = axes('Parent', figure1);
N_phase = 400;  % Grid resolution

% Define coordinate ranges
x = -2:4/(N_phase-1):2;
y = -2:4/(N_phase-1):2;

% Initialize phase arrays
phase_sign = zeros(N_phase, N_phase);
phase_sign0 = zeros(N_phase, N_phase);
tpp = zeros(N_phase, N_phase);

% Calculate phase patterns
for i = 1:N_phase
    for j = 1:N_phase
        phase_sign(i,j) = sign(cos((x(1,i)+y(1,j))*pi/4) * cos((x(1,i)-y(1,j))*pi/4));
        phase_sign0(i,j) = sign(sin((x(1,i)+y(1,j))*pi/4) * sin((x(1,i)-y(1,j))*pi/4));
    end
end

% Combine patterns to create phase regions
for i = 1:N_phase
    for j = 1:N_phase
        if phase_sign(i,j) == 1 && phase_sign0(i,j) == 1
            tpp(i,j) = 1;
        elseif phase_sign(i,j) == 1 && phase_sign0(i,j) == -1
            tpp(i,j) = 2;
        elseif phase_sign(i,j) == -1 && phase_sign0(i,j) == 1
            tpp(i,j) = 3;
        elseif phase_sign(i,j) == -1 && phase_sign0(i,j) == -1
            tpp(i,j) = 4;
        end
    end
end

% Define custom colormap
cmap = [0.5, 0.5, 0.5];  % Gray color
custom_colormap = [194/256, 206/256, 220/256;   % Light blue
                   145/256, 173/256, 158/256;   % Light green
                   216/256, 156/256, 122/256;   % Orange
                   0.5, 0, 0.5];               % Purple

% Plot phase surface
colormap(custom_colormap);
caxis([1, 4]);
hold on
surf(x, y, tpp, 'FaceColor', cmap);
xlabel('\delta/\Omega');
ylabel('t_B/\Omega_0');
zlabel('T');
shading interp;  % Remove grid lines, smooth colors

% Part 2: Parameter setting and initialization
l = 10;
grid_size = 16;        % Grid dimension
m = grid_size * grid_size;   % Number of samples
A = zeros(2, m);       % Coupling coefficients
N = 40;                % Number of momentum points
omg = 20;              % Number of quasi-energy levels
T = 2*pi/4.4;          % Driving period
freq = 2*pi/T;         % Frequency
g = 0.5*pi/T;          % Coupling strength
epsilon = 0.03;        % Variance parameter (Gaussian kernel width)

% Set up parameter grid
for i = 1:grid_size
    for j = 1:grid_size
        A(1, (i-1)*grid_size+j) = 4*pi*((j)/(grid_size+1)-0.5)/T;
        A(2, (i-1)*grid_size+j) = 4*pi*((i)/(grid_size+1.5)-0.5)/T;
    end
end

% Preallocate arrays
psi_1 = zeros(N, 2*(2*omg+1), m);
psi_2 = zeros(N, 2*(2*omg+1), m);
psi_1_t_up = zeros(N, N, m);
psi_1_t_down = zeros(N, N, m);
psi_2_t_up = zeros(N, N, m);
psi_2_t_down = zeros(N, N, m);
Heff = zeros(2*(2*omg+1), 2*(2*omg+1));
proj0 = zeros(2, 2, N, N, m);

% Pauli matrices
sigma_x = [0, 1; 1, 0];
sigma_y = [0, -1i; 1i, 0];
sigma_z = [1, 0; 0, -1];

% Main calculation loop
for i = 1:m
    for k = 1:N
        k0 = k*2*pi/N;
        
        % Define time-dependent Hamiltonians
        H0 = 1*(-A(1,i)*((sin(k0)*sigma_x)-cos(k0)*sigma_y)+2*g*sin(k0)*sigma_z+A(2,i)*sigma_y)/2;
        H1 = 1*(-A(1,i)*((sin(k0)*sigma_x)-cos(k0)*sigma_y)-A(2,i)*sigma_y)/1i/pi;
        
        Hfu1 = -H1;
        H3 = H1/3;
        Hfu3 = -H3;
        H5 = H1/5;
        Hfu5 = -H5;
        H7 = H1/7;
        Hfu7 = -H7;
        
        % Construct extended Hamiltonian
        for j = 1:2*omg+1
            Heff(2*j-1,2*j-1) = H0(1,1) + (j)*2*pi/T + 2*omg;
            Heff(2*j,2*j) = H0(2,2) + (j)*2*pi/T + 2*omg;
            Heff(2*j-1,2*j) = H0(1,2);
            Heff(2*j,2*j-1) = H0(2,1);
        end
        
        % First harmonic coupling
        for j = 1:2*omg
            Heff(2*j-1,2*j+1) = Hfu1(1,1);
            Heff(2*j,2*j+2) = Hfu1(2,2);
            Heff(2*j-1,2*j+2) = Hfu1(1,2);
            Heff(2*j,2*j+1) = Hfu1(2,1);
            Heff(2*j+1,2*j-1) = H1(1,1);
            Heff(2*j+2,2*j) = H1(2,2);
            Heff(2*j+2,2*j-1) = H1(2,1);
            Heff(2*j+1,2*j) = H1(1,2);
        end
        
        % Third harmonic coupling
        for j = 1:2*omg-2
            Heff(2*j-1,2*j+5) = Hfu3(1,1);
            Heff(2*j,2*j+6) = Hfu3(2,2);
            Heff(2*j-1,2*j+6) = Hfu3(1,2);
            Heff(2*j,2*j+5) = Hfu3(2,1);
            Heff(2*j+5,2*j-1) = H3(1,1);
            Heff(2*j+6,2*j) = H3(2,2);
            Heff(2*j+6,2*j-1) = H3(2,1);
            Heff(2*j+5,2*j) = H3(1,2);
        end
        
        % Fifth harmonic coupling
        for j = 1:2*omg-4
            Heff(2*j-1,2*j+9) = Hfu5(1,1);
            Heff(2*j,2*j+10) = Hfu5(2,2);
            Heff(2*j-1,2*j+10) = Hfu5(1,2);
            Heff(2*j,2*j+9) = Hfu5(2,1);
            Heff(2*j+9,2*j-1) = H5(1,1);
            Heff(2*j+10,2*j) = H5(2,2);
            Heff(2*j+10,2*j-1) = H5(2,1);
            Heff(2*j+9,2*j) = H5(1,2);
        end
        
        % Seventh harmonic coupling
        for j = 1:2*omg-6
            Heff(2*j-1,2*j+13) = Hfu7(1,1);
            Heff(2*j,2*j+14) = Hfu7(2,2);
            Heff(2*j-1,2*j+14) = Hfu7(1,2);
            Heff(2*j,2*j+13) = Hfu7(2,1);
            Heff(2*j+13,2*j-1) = H7(1,1);
            Heff(2*j+14,2*j) = H7(2,2);
            Heff(2*j+14,2*j-1) = H7(2,1);
            Heff(2*j+13,2*j) = H7(1,2);
        end
        
        % Diagonalize extended Hamiltonian
        [V, E] = eigs(Heff, 2*(2*omg+1));
        e1 = diag(real(E) - (2*omg+(omg+1)*2*pi/T)*eye(2*(2*omg+1)));
        
        for j = 1:2*(2*omg+1)
            psi_1(k,j,i) = V(j, 2*omg+1);
            psi_2(k,j,i) = V(j, 2*omg+2);
        end
    end
end

% Time evolution via Fourier synthesis
for i = 1:m
    for k = 1:N
        for t = 1:N
            for j = 1:2*omg+1
                phase = exp(1i*(j-omg-1)*(freq)*t/N*2*pi/freq);
                psi_1_t_up(k,t,i) = psi_1_t_up(k,t,i) + phase * psi_1(k,2*j-1,i);
                psi_1_t_down(k,t,i) = psi_1_t_down(k,t,i) + phase * psi_1(k,2*j,i);
                psi_2_t_up(k,t,i) = psi_2_t_up(k,t,i) + phase * psi_2(k,2*j-1,i);
                psi_2_t_down(k,t,i) = psi_2_t_down(k,t,i) + phase * psi_2(k,2*j,i);
            end
        end
    end
end

% Calculate projection operators
for i = 1:m
    for kx = 1:N
        for t = 1:N
            psi1 = [psi_1_t_up(kx,t,i); psi_1_t_down(kx,t,i)];
            psi2 = [psi_2_t_up(kx,t,i); psi_2_t_down(kx,t,i)];
            proj0(:,:,kx,t,i) = 1 * (psi1 * psi1') - 1 * (psi2 * psi2');
        end
    end
end

% Calculate similarity matrix
similarity_mat = ones(m, m);
for j = 1:m
    for o = 1:m
        for kx = 1:N
            for t = 1:N
                det_val = det(proj0(:,:,kx,t,o) + proj0(:,:,kx,t,j));
                similarity_mat(o,j) = (1 - exp(-abs(det_val)^2/(epsilon^2))) * similarity_mat(o,j);
            end
        end
    end
end

% Normalize similarity matrix
z = zeros(1, m);
for i = 1:m
    z(1,i) = sum(similarity_mat(i,:));
end

P = zeros(m, m);
for i = 1:m
    for j = 1:m
        P(i,j) = similarity_mat(i,j) / sqrt(z(1,i)) / sqrt(z(1,j));
    end
end

% Eigenvalue decomposition
[V, E] = eigs(P, m);
e = diag(E);

% Hierarchical clustering
data_reduced = zeros(m, 8);
for i = 1:8
    data_reduced(:,i) = V(:,i);
end
data_cluster = data_reduced;
ZZ = linkage(data_cluster, 'single');
TT = cluster(ZZ, 'maxclust', 8);

% Create cluster map
f = zeros(grid_size, grid_size);
for i = 1:grid_size
    for j = 1:grid_size
        for cluster_num = 1:25
            if TT((i-1)*grid_size+j,1) == cluster_num
                f(i,j) = cluster_num;
            end
        end
    end
end

% Prepare normalized coordinates
B = zeros(2, m);
for i = 1:grid_size*grid_size
    B(1,i) = A(1,i)/pi*T;
    B(2,i) = A(2,i)/pi*T;
end

% Plot clusters with different colors and markers
for i = 1:grid_size
    for j = 1:grid_size
        switch f(i,j)
            case 1
                plot3(B(1,(i-1)*grid_size+j), B(2,(i-1)*grid_size+j), 10, '.', 'color', [0 1 0], 'MarkerSize', 40);
            case 2
                plot3(B(1,(i-1)*grid_size+j), B(2,(i-1)*grid_size+j), 10, '.', 'color', [0 0 1], 'MarkerSize', 40);
            case 3
                plot3(B(1,(i-1)*grid_size+j), B(2,(i-1)*grid_size+j), 10, '.', 'color', [1 0 1], 'MarkerSize', 40);
            case 4
                plot3(B(1,(i-1)*grid_size+j), B(2,(i-1)*grid_size+j), 10, '^-', 'MarkerFaceColor', [0 1 0], ...
                    'MarkerEdgeColor', [0 1 0], 'MarkerSize', 10);
            case 5
                plot3(B(1,(i-1)*grid_size+j), B(2,(i-1)*grid_size+j), 10, '^-', 'color', [0 0 1], ...
                    'MarkerFaceColor', [0 0 1], 'MarkerSize', 10);
            case 6
                plot3(B(1,(i-1)*grid_size+j), B(2,(i-1)*grid_size+j), 10, '^-', 'color', [1 0 0], ...
                    'MarkerFaceColor', [1 0 0], 'MarkerSize', 10);
            case 7
                plot3(B(1,(i-1)*grid_size+j), B(2,(i-1)*grid_size+j), 10, '.', 'color', [1 0 0], 'MarkerSize', 40);
            case 8
                plot3(B(1,(i-1)*grid_size+j), B(2,(i-1)*grid_size+j), 10, '^-', 'color', [1 0 1], ...
                    'MarkerFaceColor', [1 0 1], 'MarkerSize', 10);
        end
        hold on
    end
end

% Axis formatting
hold on
set(gca, 'FontName', 'Times New Roman', 'fontsize', 32)
axis square
set(gca, 'LineWidth', 5)
box on
set(gca, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top')
set(axes1, 'TickLabelInterpreter', 'latex', 'XTick', [-2 -1 0 1 2], 'YTick', [-2 -1 0 1 2], 'FontSize', 32)
ax = gca;
ax.XAxis.TickLabelColor = 'k';
ax.YAxis.TickLabelColor = 'k';
xlabel('$J_1T/\pi$', 'Interpreter', 'latex', 'Color', 'k')
ylabel('$J_2T/\pi$', 'Interpreter', 'latex', 'Color', 'k')

% PCA analysis
X = V(:, 1:8);
[coeff, score, latent, tsquared, explained, mu] = pca(X, 'Algorithm', 'svd');
figure;
plot(score(:,1), score(:,2), 'r*', 'MarkerSize', 10);

% Save data
%save('E:\mnist\D.mat', 'e', 'score', 'tpp', 'B', 'f', 'TT');