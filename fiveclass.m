% Unsupervised clustering of Floquet topological phases
% Model: three-step piecewise-constant Hamiltonian
%   H1 = 3 * gamma1 * sigma_x
%   H2 = 3 * gamma2 * sigma_y
%   H3 = 3 * gamma3 * sigma_z
% Each step lasts T/3 within one driving period T.
%
% This script scans J2 while keeping J1 = 0.5*pi and J3 = 0.2*pi fixed.
% Floquet eigenstates are obtained from the one-period evolution operator.
% A kernel matrix is then constructed from the FFO and visualized using
% diffusion-map coordinates followed by PCA.

clear;
close all;

%% ===================== Model and numerical parameters =====================
T = 1;                      % Driving period
J1 = 0.5 * pi;              % Fixed hopping/coupling strength in the x step
J3 = 0.2 * pi;              % Fixed hopping/coupling strength in the z step
M  = 1;                     % On-site mass parameter

Nk = 100;                   % Number of k-space grid points along each direction
Nt = 100;                   % Number of time grid points
m  = 25;                    % Number of samples with different J2 values

%% ===================== J2 sampling =====================
J2_min = 0;
J2_max = 4.95 * pi;

% The original script avoids the following selected points.
% Keep this setting unchanged to preserve the original sampling strategy.
avoid_vals = (0:5) * pi;
offset = 0.05;              % Minimum distance from the avoided points

% Uniformly sample J2 within the specified interval.
J2_list_raw = linspace(J2_min, J2_max, m);

% Move samples away from the avoided points if they are too close.
for idx = 1:m
    val = J2_list_raw(idx);

    for av = avoid_vals
        if abs(val - av) < offset
            if val >= av
                J2_list_raw(idx) = av + offset;
            else
                J2_list_raw(idx) = av - offset;
            end
        end
    end
end

J2_list = sort(J2_list_raw);
fprintf('Generated %d J2 samples in the range [%.3f, %.3f].\n', ...
    m, min(J2_list), max(J2_list));

%% ===================== Brillouin-zone and time grids =====================
% Discretize the Brillouin zone with kx, ky in [-pi, pi).
kx = linspace(-pi, pi, Nk + 1);
kx = kx(1:end-1);

ky = linspace(-pi, pi, Nk + 1);
ky = ky(1:end-1);

[KX, KY] = meshgrid(kx, ky);
KX = KX(:);
KY = KY(:);
num_k = numel(KX);

% Time grid from 0 to T.
t_vals = linspace(0, T, Nt);

%% ===================== Pauli matrices =====================
sx = [0, 1; 1, 0];
sy = [0, -1i; 1i, 0];
sz = [1, 0; 0, -1];
I2 = eye(2);

% Store wave functions as psi_all(k_index, time_index, spinor_component, sample_index).
psi_all = zeros(num_k, Nt, 2, m);

%% ===================== Floquet eigenstates and time evolution =====================
fprintf('Computing Floquet eigenstates and time-evolved wave functions...\n');

for i_sample = 1:m
    J2_val = J2_list(i_sample);

    fprintf('  Sample %d/%d: J2/pi = %.6f\n', ...
        i_sample, m, J2_val / pi);

    % Evaluate gamma functions over all k points.
    gamma1 = J1 * sin(KX);
    gamma2 = J2_val * sin(KY);
    gamma3 = J3 * (M + cos(KX) + cos(KY));

    % Use parfor for acceleration. Replace parfor with for if the Parallel
    % Computing Toolbox is unavailable.
    parfor ik = 1:num_k
        g1 = 3 * gamma1(ik);
        g2 = 3 * gamma2(ik);
        g3 = 3 * gamma3(ik);

        % Full evolution operators for the three steps, each with duration T/3.
        U1_full = expm(-1i * g1 * T * sx / 3);
        U2_full = expm(-1i * g2 * T * sy / 3);
        U3_full = expm(-1i * g3 * T * sz / 3);

        % One-period evolution operator: U(T) = U3 * U2 * U1.
        U_total = U3_full * U2_full * U1_full;

        % Diagonalize U(T) to obtain quasienergies and Floquet eigenstates.
        [V, D] = eig(U_total);
        quasi_e = -angle(diag(D)) / T;

        % Select the band with the lower quasienergy as the occupied state.
        [~, idx_occ] = min(quasi_e);
        psi0 = V(:, idx_occ);
        psi0 = psi0 / norm(psi0);

        % Compute the time-evolved wave function psi(k, t).
        psi_t = zeros(Nt, 2);

        for it = 1:Nt
            t_now = t_vals(it);

            if t_now < T / 3
                theta = g1 * t_now;
                U = cos(theta) * I2 - 1i * sin(theta) * sx;

            elseif t_now < 2 * T / 3
                theta2 = g2 * (t_now - T / 3);
                U = (cos(theta2) * I2 - 1i * sin(theta2) * sy) * U1_full;

            else
                theta3 = g3 * (t_now - 2 * T / 3);
                U = (cos(theta3) * I2 - 1i * sin(theta3) * sz) ...
                    * U2_full * U1_full;
            end

            psi_t(it, :) = (U * psi0).';
        end

        psi_all(ik, :, :, i_sample) = psi_t;
    end
end

%% ===================== Kernel matrix from FFO =====================
fprintf('Computing the continuous FFO-based kernel matrix...\n');

epsilon = 1e-3;             % Scale factor controlling the exponential decay
K = ones(m, m);             % Initialize the kernel matrix

for i = 1:m
    fprintf('  Kernel row %d/%d\n', i, m);

    for j = i+1:m
        prod_val = 1.0;

        for ik = 1:num_k
            for it = 1:Nt
                psi_i = squeeze(psi_all(ik, it, :, i));
                psi_j = squeeze(psi_all(ik, it, :, j));

                psi_i = psi_i(:);
                psi_j = psi_j(:);

                Q_i = I2 - 2 * (psi_i * psi_i');
                Q_j = I2 - 2 * (psi_j * psi_j');

                det_abs = abs(det(Q_i + Q_j));

                % Exponential kernel factor: 1 - exp(-det_abs^2 / epsilon^2).
                fac = 1 - exp(-det_abs^2 / epsilon^2);
                prod_val = prod_val * fac;
            end
        end

        K(i, j) = prod_val;
        K(j, i) = prod_val;
    end
end

%% ===================== Diffusion map =====================
fprintf('Applying diffusion-map normalization and PCA visualization...\n');

% Compute row sums of the kernel matrix.
z = sum(K, 2);

% Prevent division by zero for isolated samples.
z(z == 0) = 1;

% Symmetric normalization.
P_sym = zeros(m, m);
for i = 1:m
    for j = 1:m
        P_sym(i, j) = K(i, j) / sqrt(z(i) * z(j));
    end
end

% Eigenvalue decomposition of the normalized kernel matrix.
[V_sym, E_sym] = eig(P_sym);
[eigvals, idx] = sort(diag(E_sym), 'descend');
V_sym = V_sym(:, idx);

% Keep the original setting: use the leading eigenvectors as diffusion coordinates.
if m >= 4
    coords = V_sym(:, 1:5);
else
    coords = V_sym(:, 2:end);
end

% Project the diffusion coordinates to two dimensions using PCA.
[~, score] = pca(coords);

%% ===================== Visualization =====================
figure('Color', 'w', 'Position', [100, 100, 600, 500]);
scatter(score(:, 1), score(:, 2), 80, J2_list / pi, 'filled');

colormap(jet);
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$J_2/\pi$';
cb.Label.FontSize = 14;

xlabel('PC 1', 'FontSize', 14);
ylabel('PC 2', 'FontSize', 14);
title('Diffusion Map Clustering of Floquet Topological Phases (Scan $J_2$)', ...
    'Interpreter', 'latex', 'FontSize', 14);

set(gca, 'FontSize', 12);
grid on;
box on;

figure('Color', 'w');
imagesc(K);
colormap([1, 1, 1; 0, 0, 0]);
axis square;

title('Kernel Matrix $K_{ij}$ (Scan $J_2$)', ...
    'Interpreter', 'latex', 'FontSize', 14);
xlabel('Sample Index $i$', 'FontSize', 12);
ylabel('Sample Index $j$', 'FontSize', 12);
set(gca, 'FontSize', 12);

%% ===================== Save results =====================
result_dir = fullfile(pwd, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

save(fullfile(result_dir, '2DD1.mat'), 'eigvals', 'score', 'V_sym');

fprintf('Program finished. Results saved to %s\n', ...
    fullfile(result_dir, '2DD1.mat'));
