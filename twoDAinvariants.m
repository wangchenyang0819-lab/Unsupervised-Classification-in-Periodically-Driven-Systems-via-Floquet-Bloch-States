clear
% Physical parameters
omega = 2 * pi;          % Driving frequency
N_k = 60;                % Momentum grid points
N_t = 60;                % Time steps
T = 2 * pi / omega;      % Driving period
delta_ab = 0.5 * pi / T; % Energy offset

% Simulation parameters
n_points = 30;           % Number of parameter points
topological_invariants = zeros(2, n_points);

% Pauli matrices
sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0];
sigma_z = [1 0; 0 -1];

% Main simulation loop
parfor param_idx = 1:n_points
    % Parameter value (0 to π)
    f = 3 * param_idx / n_points;
    J = f * pi / T;
    
    % Initialize Hamiltonian components
    H0 = zeros(2, 2, N_k+1, N_k+1);
    H1 = zeros(2, 2, N_k+1, N_k+1);
    H2 = zeros(2, 2, N_k+1, N_k+1);
    H3 = zeros(2, 2, N_k+1, N_k+1);
    H_offset = zeros(2, 2, N_k+1, N_k+1);
    
    % Evolution operator arrays
    U = zeros(2, 2, N_k+1, N_k+1, N_t+1);
    
    % Construct momentum-space Hamiltonians
    for kx_idx = 1:N_k+1
        for ky_idx = 1:N_k+1
            kx = 2 * pi * (kx_idx-1) / N_k;
            ky = 2 * pi * (ky_idx-1) / N_k;
            
            % Constant term
            H0(:, :, kx_idx, ky_idx) = J * sigma_x;
            
            % Nearest-neighbor hopping terms
            H1(1, 2, kx_idx, ky_idx) = J * exp(1i * kx);
            H1(2, 1, kx_idx, ky_idx) = J * exp(-1i * kx);
            
            H2(1, 2, kx_idx, ky_idx) = J * exp(1i * (kx + ky));
            H2(2, 1, kx_idx, ky_idx) = J * exp(-1i * (kx + ky));
            
            H3(1, 2, kx_idx, ky_idx) = J * exp(1i * ky);
            H3(2, 1, kx_idx, ky_idx) = J * exp(-1i * ky);
            
            % Energy offset term
            H_offset(:, :, kx_idx, ky_idx) = delta_ab * sigma_z;
        end
    end
    
    % Time evolution (piecewise constant driving)
    for t_idx = 1:N_t/5
        for kx_idx = 1:N_k+1
            for ky_idx = 1:N_k+1
                H = H0(:, :, kx_idx, ky_idx) + H_offset(:, :, kx_idx, ky_idx);
                dt = (t_idx-1)/N_t * T;
                U(:, :, kx_idx, ky_idx, t_idx) = expm(-1i * H * dt);
            end
        end
    end
    
    for t_idx = N_t/5+1:2*N_t/5
        for kx_idx = 1:N_k+1
            for ky_idx = 1:N_k+1
                H = H1(:, :, kx_idx, ky_idx) + H_offset(:, :, kx_idx, ky_idx);
                dt = (t_idx - N_t/5)/N_t * T;
                U(:, :, kx_idx, ky_idx, t_idx) = expm(-1i * H * dt) * U(:, :, kx_idx, ky_idx, N_t/5);
            end
        end
    end
    
    for t_idx = 2*N_t/5+1:3*N_t/5
        for kx_idx = 1:N_k+1
            for ky_idx = 1:N_k+1
                H = H2(:, :, kx_idx, ky_idx) + H_offset(:, :, kx_idx, ky_idx);
                dt = (t_idx - 2*N_t/5)/N_t * T;
                U(:, :, kx_idx, ky_idx, t_idx) = expm(-1i * H * dt) * U(:, :, kx_idx, ky_idx, 2*N_t/5);
            end
        end
    end
    
    for t_idx = 3*N_t/5+1:4*N_t/5
        for kx_idx = 1:N_k+1
            for ky_idx = 1:N_k+1
                H = H3(:, :, kx_idx, ky_idx) + H_offset(:, :, kx_idx, ky_idx);
                dt = (t_idx - 3*N_t/5)/N_t * T;
                U(:, :, kx_idx, ky_idx, t_idx) = expm(-1i * H * dt) * U(:, :, kx_idx, ky_idx, 3*N_t/5);
            end
        end
    end
    
    for t_idx = 4*N_t/5+1:N_t
        for kx_idx = 1:N_k+1
            for ky_idx = 1:N_k+1
                H = H_offset(:, :, kx_idx, ky_idx);
                dt = (t_idx - 4*N_t/5)/N_t * T;
                U(:, :, kx_idx, ky_idx, t_idx) = expm(-1i * H * dt) * U(:, :, kx_idx, ky_idx, 4*N_t/5);
            end
        end
    end
    
    % Diagonalize final time evolution operator
    eigenvalues = zeros(2, N_k+1, N_k+1);
    eigenvectors = zeros(2, 2, N_k+1, N_k+1);
    
    for kx_idx = 1:N_k+1
        for ky_idx = 1:N_k+1
            [V, D] = eig(U(:, :, kx_idx, ky_idx, N_t));
            eigenvalues(:, kx_idx, ky_idx) = diag(D);
            eigenvectors(:, :, kx_idx, ky_idx) = V;
        end
    end
    
    % Compute effective Hamiltonians in 0 and π gap
    angle_pi = zeros(2, N_k+1, N_k+1);
    angle_0 = zeros(2, N_k+1, N_k+1);
    H_eff_0 = zeros(2, 2, N_k+1, N_k+1);
    H_eff_pi = zeros(2, 2, N_k+1, N_k+1);
    
    for kx_idx = 1:N_k+1
        for ky_idx = 1:N_k+1
            % Quasi-energies (shifted by 2π/T for π gap)
            angle_pi(1, kx_idx, ky_idx) = -angle(eigenvalues(1, kx_idx, ky_idx)) / T - 2*pi/T;
            angle_pi(2, kx_idx, ky_idx) = -angle(eigenvalues(2, kx_idx, ky_idx)) / T - 2*pi/T;
            
            % Map to principal branch (-π/T, π/T)
            angle_0(1, kx_idx, ky_idx) = mod(angle_pi(1, kx_idx, ky_idx) + 2*pi/T, 2*pi/T) - pi/T;
            angle_0(2, kx_idx, ky_idx) = mod(angle_pi(2, kx_idx, ky_idx) + 2*pi/T, 2*pi/T) - pi/T;
            
            % Construct effective Hamiltonians
            H_eff_0(:, :, kx_idx, ky_idx) = angle_0(1, kx_idx, ky_idx) * (eigenvectors(:, 1, kx_idx, ky_idx) * eigenvectors(:, 1, kx_idx, ky_idx)') + ...
                                           angle_0(2, kx_idx, ky_idx) * (eigenvectors(:, 2, kx_idx, ky_idx) * eigenvectors(:, 2, kx_idx, ky_idx)');
            H_eff_pi(:, :, kx_idx, ky_idx) = angle_pi(1, kx_idx, ky_idx) * (eigenvectors(:, 1, kx_idx, ky_idx) * eigenvectors(:, 1, kx_idx, ky_idx)') + ...
                                            angle_pi(2, kx_idx, ky_idx) * (eigenvectors(:, 2, kx_idx, ky_idx) * eigenvectors(:, 2, kx_idx, ky_idx)');
        end
    end
    
    % Construct modified evolution operators
    U_modified_0 = zeros(2, 2, N_k+1, N_k+1, N_t+1);
    U_modified_pi = zeros(2, 2, N_k+1, N_k+1, N_t+1);
    
    for t_idx = 1:N_t
        for kx_idx = 1:N_k+1
            for ky_idx = 1:N_k+1
                dt = (t_idx-1)/N_t * T;
                U_modified_0(:, :, kx_idx, ky_idx, t_idx) = U(:, :, kx_idx, ky_idx, t_idx) * expm(1i * H_eff_0(:, :, kx_idx, ky_idx) * dt);
                U_modified_pi(:, :, kx_idx, ky_idx, t_idx) = U(:, :, kx_idx, ky_idx, t_idx) * expm(1i * H_eff_pi(:, :, kx_idx, ky_idx) * dt);
            end
        end
    end
    
    % Calculate topological invariants (3D winding numbers)
    winding_0 = 0;
    winding_pi = 0;
    
    for t_idx = 1:N_t-1
        for kx_idx = 1:N_k
            for ky_idx = 1:N_k
                % For 0 gap
                Ut = U_modified_0(:, :, kx_idx, ky_idx, t_idx);
                Ut_inv = inv(Ut);
                
                dUdt = Ut_inv * (U_modified_0(:, :, kx_idx, ky_idx, t_idx+1) - U_modified_0(:, :, kx_idx, ky_idx, t_idx));
                dUdx = Ut_inv * (U_modified_0(:, :, kx_idx+1, ky_idx, t_idx) - U_modified_0(:, :, kx_idx, ky_idx, t_idx));
                dUdy = Ut_inv * (U_modified_0(:, :, kx_idx, ky_idx+1, t_idx) - U_modified_0(:, :, kx_idx, ky_idx, t_idx));
                
                comm = dUdx*dUdy*dUdt - dUdx*dUdt*dUdy - dUdy*dUdx*dUdt + ...
                       dUdy*dUdt*dUdx + dUdt*dUdx*dUdy - dUdt*dUdy*dUdx;
                winding_0 = winding_0 + trace(comm) / (24 * pi^2);
                
                % For π gap
                Ut_pi = U_modified_pi(:, :, kx_idx, ky_idx, t_idx);
                Ut_pi_inv = inv(Ut_pi);
                
                dUdt_pi = Ut_pi_inv * (U_modified_pi(:, :, kx_idx, ky_idx, t_idx+1) - U_modified_pi(:, :, kx_idx, ky_idx, t_idx));
                dUdx_pi = Ut_pi_inv * (U_modified_pi(:, :, kx_idx+1, ky_idx, t_idx) - U_modified_pi(:, :, kx_idx, ky_idx, t_idx));
                dUdy_pi = Ut_pi_inv * (U_modified_pi(:, :, kx_idx, ky_idx+1, t_idx) - U_modified_pi(:, :, kx_idx, ky_idx, t_idx));
                
                comm_pi = dUdx_pi*dUdy_pi*dUdt_pi - dUdx_pi*dUdt_pi*dUdy_pi - dUdy_pi*dUdx_pi*dUdt_pi + ...
                          dUdy_pi*dUdt_pi*dUdx_pi + dUdt_pi*dUdx_pi*dUdy_pi - dUdt_pi*dUdy_pi*dUdx_pi;
                winding_pi = winding_pi + trace(comm_pi) / (24 * pi^2);
            end
        end
    end
    
    % Store results
    topological_invariants(:, param_idx) = [real(winding_0); real(winding_pi)];
end

% Plot results
figure('Position', [100, 100, 800, 600]);
imagesc(topological_invariants);
colorbar;
xlabel('Parameter Index (aa)');
ylabel('Topological Gap (0/\pi)');
title('Topological Invariants vs Parameter');
set(gca, 'YTick', [1, 2], 'YTickLabel', {'0-gap', '\pi-gap'});
colormap(jet);