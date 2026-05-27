% Unsupervised clustering of Floquet topological phases.
%
% Model:
%   H1 = 3 * gamma1 * sigma_x
%   H2 = 3 * gamma2 * sigma_y
%   H3 = 3 * gamma3 * sigma_z
%
% Each segment lasts T/3 within one driving period T.
% Floquet eigenstates are obtained by diagonalizing the time-evolution
% operator. The FFO-based kernel matrix is then used for diffusion-map
% clustering.

clear;
close all;

%% Parameter settings

T = 1;                       % Driving period.
J1 = pi / 3;                 % Hopping parameter along the x direction.
J2 = pi / 3;                 % Hopping parameter along the y direction.
M = 1;                       % On-site mass parameter.

Nk = 100;                    % Number of k-space grid points in each direction.
Nt = 100;                    % Number of time grid points.
m = 28;                      % Number of samples with different J3 values.

%% Uniformly sample J3 while avoiding transition points

J3_min = 0;
J3_max = 7 * pi / 3;
avoid_vals = (0:6.95) * pi / 3;   % Transition points.
offset = 0.05;                    % Exclusion width around transition points.

J3_list_raw = linspace(J3_min + 0.01, J3_max - 0.01, m);

for idx = 1:m
    val = J3_list_raw(idx);

    for av = avoid_vals
        if abs(val - av) < offset
            if val >= av
                J3_list_raw(idx) = av + offset;
            else
                J3_list_raw(idx) = av - offset;
            end
        end
    end
end

J3_list = sort(J3_list_raw);

%% Discretize the Brillouin zone

kx = linspace(-pi, pi, Nk + 1);
kx = kx(1:end - 1);

ky = linspace(-pi, pi, Nk + 1);
ky = ky(1:end - 1);

[KX, KY] = meshgrid(kx, ky);
KX = KX(:);
KY = KY(:);

num_k = length(KX);

%% Discretize one driving period

t_vals = linspace(0, T, Nt);

%% Pauli matrices

sx = [0, 1; 1, 0];
sy = [0, -1i; 1i, 0];
sz = [1, 0; 0, -1];
I2 = eye(2);

%% Compute Floquet eigenstates and time-evolved wavefunctions

% psi_all stores psi(kx, ky, t, spinor_component, sample_index).
psi_all = zeros(num_k, Nt, 2, m);

fprintf('Computing Floquet eigenstates and time-evolved wavefunctions...\n');

for i_sample = 1:m
    J3_val = J3_list(i_sample);

    gamma1 = J1 * sin(KX);
    gamma2 = J2 * sin(KY);
    gamma3 = J3_val * (M + cos(KX) + cos(KY));

    % Use a regular "for" loop instead of "parfor" if the Parallel
    % Computing Toolbox is unavailable.
    parfor ik = 1:num_k
        g1 = 3 * gamma1(ik);
        g2 = 3 * gamma2(ik);
        g3 = 3 * gamma3(ik);

        % Full evolution matrices for the three T/3 segments.
        U1_full = expm(-1i * g1 * T * sx / 3);
        U2_full = expm(-1i * g2 * T * sy / 3);
        U3_full = expm(-1i * g3 * T * sz / 3);

        % Full-period time-evolution operator: U(T) = U3 * U2 * U1.
        U_total = U3_full * U2_full * U1_full;

        % Diagonalize U(T) to obtain quasienergies and eigenstates.
        [V, D] = eig(U_total);
        quasi_e = -angle(diag(D)) / T;

        % Select the occupied state as the band with the lower quasienergy.
        [~, idx_occ] = min(quasi_e);
        psi0 = V(:, idx_occ);
        psi0 = psi0 / norm(psi0);

        % Compute psi(k, t) for the selected Floquet state.
        psi_t = zeros(Nt, 2);

        for it = 1:Nt
            t_now = t_vals(it);

            if t_now < T / 3
                theta = 3 * g1 * t_now;
                U = cos(theta) * I2 - 1i * sin(theta) * sx;

            elseif t_now < 2 * T / 3
                theta2 = 3 * g2 * (t_now - T / 3);
                U = (cos(theta2) * I2 - 1i * sin(theta2) * sy) * U1_full;

            else
                theta3 = 3 * g3 * (t_now - 2 * T / 3);
                U = (cos(theta3) * I2 - 1i * sin(theta3) * sz) * U2_full * U1_full;
            end

            psi_t(it, :) = (U * psi0).';
        end

        psi_all(ik, :, :, i_sample) = psi_t;
    end
end

%% Compute the FFO-based kernel matrix

fprintf('Computing the continuous FFO kernel matrix...\n');

epsilon = 1e-5;             % Scale factor controlling exponential decay.
K = ones(m, m);

for i = 1:m
    fprintf('Processing sample %d of %d...\n', i, m);

    for j = i + 1:m
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

                % Exponential kernel factor:
                %   1 - exp(-det_abs^2 / epsilon^2)
                fac = 1 - exp(-det_abs^2 / epsilon^2);
                prod_val = prod_val * fac;
            end
        end

        K(i, j) = prod_val;
        K(j, i) = prod_val;
    end
end

%% Diffusion-map embedding

fprintf('Running diffusion-map embedding and visualization...\n');

z = sum(K, 2);
z(z == 0) = 1;

P_sym = zeros(m, m);

for i = 1:m
    for j = 1:m
        P_sym(i, j) = K(i, j) / sqrt(z(i) * z(j));
    end
end

[V_sym, E_sym] = eig(P_sym);
[eigvals, idx] = sort(diag(E_sym), 'descend');
V_sym = V_sym(:, idx);

if m >= 4
    coords = V_sym(:, 1:min(7, m));
else
    coords = V_sym(:, 1:end);
end

%% PCA visualization of diffusion coordinates

[~, score] = pca(coords);

figure('Color', 'w', 'Position', [100, 100, 600, 500]);
scatter(score(:, 1), score(:, 2), 80, J3_list / pi, 'filled');

colormap(jet);
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$J_3/\pi$';
cb.Label.FontSize = 14;

xlabel('PC 1', 'FontSize', 14);
ylabel('PC 2', 'FontSize', 14);
title('Diffusion Map Clustering of Floquet Topological Phases', 'FontSize', 14);

set(gca, 'FontSize', 12);
grid on;
box on;

%% Plot the kernel matrix

figure('Color', 'w');
imagesc(K);

colormap([1, 1, 1; 0, 0, 0]);
axis square;

title('Kernel Matrix $K_{ij}$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Sample index $i$', 'FontSize', 12);
ylabel('Sample index $j$', 'FontSize', 12);

set(gca, 'FontSize', 12);

%% Save results

output_dir = fullfile(pwd, 'results');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

output_file = fullfile(output_dir, '2DD2.mat');
save(output_file, 'eigvals', 'score', 'V_sym');

fprintf('Program completed. Results saved to: %s\n', output_file);
