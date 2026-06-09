clear

% Parameters
l = 10;                     % Dimension for phase diagram
m = 1;                      % Number of samples
A = zeros(2, m);            % Coupling coefficients
N = 13;                     % Number of discrete grid points
omg = 30;                   % Number of quasi-energy levels
T = 2*pi/(1*pi);            % Driving period
freq = 2*pi/T;              % Frequency
g = 1.6;                    % Coupling strength

% Set coupling coefficients
for i = 1:m
    A(1, i) = 1.5*pi/T;     % Coupling coefficient 1
    A(2, i) = 1*pi/T;       % Coupling coefficient 2
end

% Initialize wave function arrays
psi_1 = zeros(N, 2*(2*omg+1), m);
psi_2 = zeros(N, 2*(2*omg+1), m);
psi_1_t_up = zeros(N, N, m);
psi_1_t_down = zeros(N, N, m);
psi_2_t_up = zeros(N, N, m);
psi_2_t_down = zeros(N, N, m);
Heff = zeros(2*(2*omg+1), 2*(2*omg+1));

% Pauli matrices
sigma_x = [0, 1; 1, 0];
sigma_y = [0, -1i; 1i, 0];
sigma_z = [1, 0; 0, -1];

% Main calculation loop
for i = 1:m
    for k = 1:N
        k0 = (k-1)*2*pi/(N-1);
        
        % Define Hamiltonians
        H0 = 1*(-A(1,i)*((sin(k0)*sigma_x) - cos(k0)*sigma_y) + 2*g*sin(k0)*sigma_z + A(2,i)*sigma_y)/2;
        H1 = 1*(-A(1,i)*((sin(k0)*sigma_x) - cos(k0)*sigma_y) - A(2,i)*sigma_y)/(1i*pi);
        
        Hfu1 = -H1;
        H3 = H1/3;
        Hfu3 = -H3;
        H5 = H1/5;
        Hfu5 = -H5;
        H7 = H1/7;
        Hfu7 = -H7;
        
        % Construct effective Hamiltonian
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
        
        % Diagonalize effective Hamiltonian
        [V, ~] = eigs(Heff, 2*(2*omg+1));
        
        for j = 1:2*(2*omg+1)
            psi_1(k, j, i) = V(j, 2*omg+2);
            psi_2(k, j, i) = V(j, 2*omg+1);
        end
    end
end

% Time evolution of wave functions
for i = 1:m
    for k = 1:N
        for t = 1:N
            for j = 1:2*omg+1
                phase = exp(1i*(j-omg-1)*(freq)*(t-1)/(N-1)*2*pi/freq);
                psi_1_t_up(k, t, i) = psi_1_t_up(k, t, i) + phase * psi_1(k, 2*j-1, i);
                psi_1_t_down(k, t, i) = psi_1_t_down(k, t, i) + phase * psi_1(k, 2*j, i);
                psi_2_t_up(k, t, i) = psi_2_t_up(k, t, i) + phase * psi_2(k, 2*j-1, i);
                psi_2_t_down(k, t, i) = psi_2_t_down(k, t, i) + phase * psi_2(k, 2*j, i);
            end
        end
    end
end

% Calculate Bloch vector components (xx, yy, zz kept as storage names)
xx = zeros(N, N);
yy = zeros(N, N);
zz = zeros(N, N);

for i = 1:m
    for k = 1:N
        for t = 1:N
            psi_vec = [psi_1_t_up(k, t, i); psi_1_t_down(k, t, i)];
            xx(k, t) = real(psi_vec' * sigma_x * psi_vec);
            yy(k, t) = real(psi_vec' * sigma_y * psi_vec);
            zz(k, t) = real(psi_vec' * sigma_z * psi_vec);
        end
    end
end

% Calculate Chern number
chern = 0;
for k = 1:(N-1)/2
    for t = 1:(N-1)
        vec1 = [xx(k, t), yy(k, t), zz(k, t)];
        vec2 = [xx(k+1, t), yy(k+1, t), zz(k+1, t)] - vec1;
        vec3 = [xx(k, t+1), yy(k, t+1), zz(k, t+1)] - vec1;
        berry_contrib = 1/(4*pi) * dot(vec1, cross(vec2, vec3)) / (2*pi/N) / (T/N);
        chern = chern + berry_contrib * (2*pi/N) * (T/N);
    end
end

% Display results
fprintf('Chern number: %f\n', chern);
fprintf('xx, yy, zz calculated. Matrix size: %d x %d\n', N, N);

% Save data
save('E:\mnist\kernelD151.mat', 'xx', 'yy', 'zz');