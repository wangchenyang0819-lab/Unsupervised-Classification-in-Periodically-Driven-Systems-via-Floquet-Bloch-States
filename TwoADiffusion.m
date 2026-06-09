clear

% Parameters
l = 10;                     % Dimension for phase diagram
m = 30;                     % Number of samples
A = zeros(2, m);            % Coupling coefficients
N = 40;                     % Number of momentum grid points
Nt = 40;                    % Number of time grid points
omg = 7;                    % Number of quasi-energy levels
T = 2*pi/(1*pi);            % Driving period
freq = 2*pi/T;              % Frequency
g = 1*pi/T;                 % Coupling strength

% Initialize wave function arrays
psi_0 = zeros(N, N, 2*(2*omg+1), m);
psi_1 = zeros(N, N, 2*(2*omg+1), m);
psi_2 = zeros(N, N, 2*(2*omg+1), m);
psi_0_t_up = zeros(N, N, Nt, m);
psi_0_t_down = zeros(N, N, Nt, m);
psi_1_t_up = zeros(N, N, Nt, m);
psi_1_t_down = zeros(N, N, Nt, m);
psi_2_t_up = zeros(N, N, Nt, m);
psi_2_t_down = zeros(N, N, Nt, m);

% Pauli matrices
sigma_x = [0, 1; 1, 0];
sigma_y = [0, -1i; 1i, 0];
sigma_z = [1, 0; 0, -1];
sigma_plus = (sigma_x + 1i*sigma_y)/2;
sigma_minus = (sigma_x - 1i*sigma_y)/2;

% Set coupling coefficients for different samples
for i = 1:m
    A(1,i) = (3*((i)/(m+1)))*pi/T;  % Hopping amplitude
    A(2,i) = 0.5*pi/T;              % Sublattice potential
end

% Main calculation loop for each sample
for i = 1:m
    for kx = 1:N
        for ky = 1:N
            kx0 = kx*2*pi/N;
            ky0 = ky*2*pi/N;
            
            % Define static Hamiltonian
            H0 = A(1,i)*(1+exp(1i*ky0)+exp(1i*(kx0+ky0))+exp(1i*kx0))*sigma_plus/5 + ...
                 A(2,i)*sigma_z + ...
                 A(1,i)*(1+exp(-1i*ky0)+exp(-1i*(kx0+ky0))+exp(-1i*kx0))*sigma_minus/5;
            
            % Compute Fourier components of driving Hamiltonian
            H1 = computeH(A(1,i), kx0, ky0, 1, sigma_plus, sigma_minus);
            Hfu1 = computeH(A(1,i), kx0, ky0, -1, sigma_plus, sigma_minus);
            H2 = computeH(A(1,i), kx0, ky0, 2, sigma_plus, sigma_minus);
            Hfu2 = computeH(A(1,i), kx0, ky0, -2, sigma_plus, sigma_minus);
            H3 = computeH(A(1,i), kx0, ky0, 3, sigma_plus, sigma_minus);
            Hfu3 = computeH(A(1,i), kx0, ky0, -3, sigma_plus, sigma_minus);
            H4 = computeH(A(1,i), kx0, ky0, 4, sigma_plus, sigma_minus);
            Hfu4 = computeH(A(1,i), kx0, ky0, -4, sigma_plus, sigma_minus);
            H5 = computeH(A(1,i), kx0, ky0, 5, sigma_plus, sigma_minus);
            Hfu5 = computeH(A(1,i), kx0, ky0, -5, sigma_plus, sigma_minus);
            
            % Build effective (Floquet) Hamiltonian
            Heff = zeros(2*(2*omg+1), 2*(2*omg+1));
            
            for j = 1:2*omg+1
                Heff(2*j-1,2*j-1) = H0(1,1) + (j)*2*pi/T + 2*omg;
                Heff(2*j,2*j) = H0(2,2) + (j)*2*pi/T + 2*omg;
                Heff(2*j-1,2*j) = H0(1,2);
                Heff(2*j,2*j-1) = H0(2,1);
            end
            
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
            
            for j = 1:2*omg-1
                Heff(2*j-1,2*j+3) = Hfu2(1,1);
                Heff(2*j,2*j+4) = Hfu2(2,2);
                Heff(2*j-1,2*j+4) = Hfu2(1,2);
                Heff(2*j,2*j+3) = Hfu2(2,1);
                Heff(2*j+3,2*j-1) = H2(1,1);
                Heff(2*j+4,2*j) = H2(2,2);
                Heff(2*j+4,2*j-1) = H2(2,1);
                Heff(2*j+3,2*j) = H2(1,2);
            end
            
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
            
            for j = 1:2*omg-3
                Heff(2*j-1,2*j+7) = Hfu4(1,1);
                Heff(2*j,2*j+8) = Hfu4(2,2);
                Heff(2*j-1,2*j+8) = Hfu4(1,2);
                Heff(2*j,2*j+7) = Hfu4(2,1);
                Heff(2*j+7,2*j-1) = H4(1,1);
                Heff(2*j+8,2*j) = H4(2,2);
                Heff(2*j+8,2*j-1) = H4(2,1);
                Heff(2*j+7,2*j) = H4(1,2);
            end
            
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
            
            % Diagonalize effective Hamiltonian
            [V, E] = eigs(Heff, 2*(2*omg+1));
            e1 = diag(real(E));
            
            for j = 1:2*(2*omg+1)
                psi_0(kx,ky,j,i) = V(j, 2*omg);
                psi_1(kx,ky,j,i) = V(j, 2*omg+1);
                psi_2(kx,ky,j,i) = V(j, 2*omg+2);
            end
        end
    end
end

% Time evolution of wave functions
for i = 1:m
    for kx = 1:N
        for ky = 1:N
            for t = 1:Nt
                for j = 1:2*omg+1
                    phase = exp(1i*(j-omg-1)*(freq)*(t-1)/(Nt-1)*2*pi/freq);
                    psi_0_t_up(kx,ky,t,i) = psi_0_t_up(kx,ky,t,i) + phase * psi_0(kx,ky,2*j-1,i);
                    psi_0_t_down(kx,ky,t,i) = psi_0_t_down(kx,ky,t,i) + phase * psi_0(kx,ky,2*j,i);
                    psi_1_t_up(kx,ky,t,i) = psi_1_t_up(kx,ky,t,i) + phase * psi_1(kx,ky,2*j-1,i);
                    psi_1_t_down(kx,ky,t,i) = psi_1_t_down(kx,ky,t,i) + phase * psi_1(kx,ky,2*j,i);
                    psi_2_t_up(kx,ky,t,i) = psi_2_t_up(kx,ky,t,i) + phase * psi_2(kx,ky,2*j-1,i);
                    psi_2_t_down(kx,ky,t,i) = psi_2_t_down(kx,ky,t,i) + phase * psi_2(kx,ky,2*j,i);
                end
            end
        end
    end
end

% Calculate projection operators
proj0 = zeros(2, 2, N, N, Nt, m);
proj_pi = zeros(2, 2, N, N, Nt, m);

for i = 1:m
    fprintf('Processing sample %d\n', i);
    for kx = 1:N
        for ky = 1:N
            for t = 1:Nt
                psi1 = [psi_1_t_up(kx,ky,t,i); psi_1_t_down(kx,ky,t,i)];
                psi2 = [psi_2_t_up(kx,ky,t,i); psi_2_t_down(kx,ky,t,i)];
                psi0 = [psi_0_t_up(kx,ky,t,i); psi_0_t_down(kx,ky,t,i)];
                
                proj0(:,:,kx,ky,t,i) = 1 * (psi1 * psi1') - 1 * (psi2 * psi2');
                proj_pi(:,:,kx,ky,t,i) = e1(2*omg+1,1) * (psi1 * psi1') + (2*pi/T - e1(2*omg+1,1)) * (psi0 * psi0');
            end
        end
    end
end

% Calculate similarity matrix between samples
similarity_mat = ones(m, m);
epsilon = 0.1;

for j = 1:m
    fprintf('Processing similarity for sample %d\n', j);
    for o = 1:m
        for kx = 1:N
            for ky = 1:N
                for t = 1:Nt
                    det_val = det(proj0(:,:,kx,ky,t,o) + proj0(:,:,kx,ky,t,j));
                    similarity_mat(o,j) = (1 - exp(-abs(det_val)^2/(epsilon^2))) * similarity_mat(o,j);
                end
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

% Plotting
theta = 0.1;
Neta = m;
etalist = linspace(0, 3, Neta);
NOmega = m;
Omegalist = linspace(0, 3, NOmega);
W = similarity_mat;

figure;
axes1 = axes;
hold(axes1, 'on');

% Define custom colormap
custom_cmap = [
    5, 48, 97;
    255, 255, 255;
    105, 0, 31] / 255;
n_colors = 256;
custom_cmap_interp = interp1(linspace(0, 1, size(custom_cmap, 1)), custom_cmap, linspace(0, 1, n_colors));
colormap(custom_cmap_interp);

% Create heatmap
for nnx = 1:Neta-1
    for nny = 1:NOmega-1
        K1a = etalist(nnx);  K1b = Omegalist(nny);
        K2a = etalist(nnx+1);  K2b = Omegalist(nny);
        K3a = etalist(nnx+1);  K3b = Omegalist(nny+1);
        K4a = etalist(nnx);  K4b = Omegalist(nny+1);
        
        fill([K1a; K2a; K3a; K4a], [K1b; K2b; K3b; K4b], ...
             [W(nnx,nny); W(nnx+1,nny); W(nnx+1,nny+1); W(nnx,nny+1)], 'edgecolor', 'none');
    end
end

box(axes1, 'on');
axis(axes1, 'square');
hold(axes1, 'off');

% Axis settings
set(axes1, 'TickLabelInterpreter', 'latex', 'XTick', [0 1 2 3], 'YTick', [0 1 2 3], 'FontSize', 32);
colorbar(axes1, 'Ticks', [0 0.5 0.999], 'TickLabels', {'$0$', '$0.5$', '$1$'}, ...
    'TickLabelInterpreter', 'latex', 'FontSize', 32, 'Color', [0 0 0]);

xlim([0 3]);
ylim([0 3]);
axis square;
set(gca, 'LineWidth', 5);
box on;
set(gca, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top');
ax = gca;
ax.XAxis.TickLabelColor = 'k';
ax.YAxis.TickLabelColor = 'k';
xlabel('$JT/\pi$', 'Interpreter', 'latex', 'fontsize', 32, 'Color', 'k');
ylabel('$JT/\pi$', 'Interpreter', 'latex', 'fontsize', 32, 'Color', 'k');

% PCA analysis
X = V(:, 1:3);
[coeff, score, latent, tsquared, explained, mu] = pca(X, 'Algorithm', 'svd');

figure;
plot(score(:,1), score(:,2), 'r*', 'MarkerSize', 10);

% Save data
%save('E:\mnist\A1.mat', 'e', 'score');

% Helper function: compute Fourier component of driving Hamiltonian
function H = computeH(A, kx0, ky0, n, sigma_plus, sigma_minus)
    phase_coeff = 1i/(2*pi*n) * A;
    term1 = (exp(-1i*2*pi/5*n) - 1) * (sigma_plus + sigma_minus);
    term2 = (exp(-1i*4*pi/5*n) - exp(-1i*2*pi/5*n)) * (exp(1i*kx0)*sigma_plus + exp(-1i*kx0)*sigma_minus);
    term3 = (exp(-1i*6*pi/5*n) - exp(-1i*4*pi/5*n)) * (exp(1i*(kx0+ky0))*sigma_plus + exp(-1i*(kx0+ky0))*sigma_minus);
    term4 = (exp(-1i*8*pi/5*n) - exp(-1i*6*pi/5*n)) * (exp(1i*ky0)*sigma_plus + exp(-1i*ky0)*sigma_minus);
    H = phase_coeff * (term1 + term2 + term3 + term4);
end