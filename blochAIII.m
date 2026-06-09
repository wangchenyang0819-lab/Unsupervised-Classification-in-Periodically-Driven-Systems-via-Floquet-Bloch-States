clear

% Parameters
N = 13;                     % Number of discrete grid points
omg = 30;                   % Number of quasi-energy levels
T = 2*pi/(1*pi);            % Driving period
freq = 2*pi/T;              % Frequency
g = 1.6;                    % Coupling strength
A22 = 0.3*pi;               % Additional coupling parameter
A1 = 0.3*pi;                % Coupling coefficient 1
A2 = 0.1*pi;                % Coupling coefficient 2

% Initialize wave function arrays
psi_1 = zeros(N, 2*(2*omg+1));
psi_1_t_up = zeros(N, N);
psi_1_t_down = zeros(N, N);

% Pauli matrices
sigma_x = [0, 1; 1, 0];
sigma_y = [0, -1i; 1i, 0];
sigma_z = [1, 0; 0, -1];

% Main calculation loop
for k = 1:N
    k0 = (k-1)*2*pi/(N-1);
    
    % Define Hamiltonians
    H0 = A2*sigma_y + (A1*sigma_x + A22*((cos(k0)*sigma_x) + (sin(k0))*sigma_y));
    H1 = g*((cos(k0)*sigma_x) + (sin(k0))*sigma_y);
    
    % Construct effective Hamiltonian
    Heff = zeros(2*(2*omg+1), 2*(2*omg+1));
    
    for j = 1:2*omg+1
        Heff(2*j-1, 2*j-1) = H0(1, 1) + (j)*2*pi/T + 2*omg;
        Heff(2*j, 2*j) = H0(2, 2) + (j)*2*pi/T + 2*omg;
        Heff(2*j-1, 2*j) = H0(1, 2);
        Heff(2*j, 2*j-1) = H0(2, 1);
    end
    
    for j = 1:2*omg
        Heff(2*j-1, 2*j+1) = H1(1, 1);
        Heff(2*j, 2*j+2) = H1(2, 2);
        Heff(2*j-1, 2*j+2) = H1(1, 2);
        Heff(2*j, 2*j+1) = H1(2, 1);
        Heff(2*j+1, 2*j-1) = H1(1, 1);
        Heff(2*j+2, 2*j) = H1(2, 2);
        Heff(2*j+2, 2*j-1) = H1(2, 1);
        Heff(2*j+1, 2*j) = H1(1, 2);
    end
    
    % Diagonalize effective Hamiltonian
    [V, ~] = eigs(Heff, 2*(2*omg+1));
    
    for j = 1:2*(2*omg+1)
        psi_1(k, j) = V(j, 2*omg+2);
    end
end

% Time evolution of wave functions
for k = 1:N
    for t = 1:N
        for j = 1:2*omg+1
            phase = exp(1i*(j-omg-1)*(freq)*(t-1)/(N-1)*2*pi/freq);
            psi_1_t_up(k, t) = psi_1_t_up(k, t) + phase * psi_1(k, 2*j-1);
            psi_1_t_down(k, t) = psi_1_t_down(k, t) + phase * psi_1(k, 2*j);
        end
    end
end

% Calculate Bloch vector components
xx = zeros(N, N);
yy = zeros(N, N);
zz = zeros(N, N);

for k = 1:N
    for t = 1:N
        psi_vec = [psi_1_t_up(k, t); psi_1_t_down(k, t)];
        xx(k, t) = real(psi_vec' * sigma_x * psi_vec);
        yy(k, t) = real(psi_vec' * sigma_y * psi_vec);
        zz(k, t) = real(psi_vec' * sigma_z * psi_vec);
    end
end
% Save data
save('E:\mnist\kernelAIIIpi.mat', 'xx', 'yy', 'zz');