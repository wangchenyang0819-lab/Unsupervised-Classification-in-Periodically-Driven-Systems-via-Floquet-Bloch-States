clear
l = 10;
m = 1;                     % Number of samples
A = zeros(2, m);           % Coupling coefficients
N = 200;                   % Number of discretization grid points
omg = 10;                  % Number of quasi-energy levels
T = 2*pi;
freq = 2*pi/T;             % Driving frequency
g = 1*pi/T;

% Start parallel pool
if isempty(gcp('nocreate'))
    parpool;               % Use default configuration to start parallel pool
end

% Preallocate large arrays
psi_1 = zeros(N, N, 2*(2*omg+1));
psi_1_t_up = zeros(N, N, N);
psi_1_t_down = zeros(N, N, N);

% Sample parameters
f = 27;%4 trival 27nontrival
A(1,1) = (3*((f)/30))*pi/T;   % Hopping amplitude
A(2,1) = 0.5*pi/T;            % Sublattice potential

% Parallelize the kx loop
parfor kx = 1:N
    % Each worker has its own copy of Pauli matrices
    sigma_x = [0, 1; 1, 0];
    sigma_y = [0, -1i; 1i, 0];
    sigma_z = [1, 0; 0, -1];
    sigma_plus = (sigma_x + 1i*sigma_y)/2;
    sigma_minus = (sigma_x - 1i*sigma_y)/2;
    
    % Preallocate temporary variables for the worker
    worker_psi_1 = zeros(N, 2*(2*omg+1));
    worker_psi_1_t_up = zeros(N, N);
    worker_psi_1_t_down = zeros(N, N);
    
    kx0 = (kx-1-N/2)*2*pi/(N-1);
    
    for ky = 1:N
        ky0 = (ky-1-N/2)*2*pi/(N-1);
        
        % Compute Hamiltonian components
        H0 = A(1,1)*(1 + exp(1i*(kx0+ky0)) + exp(1i*(2*kx0)) + exp(1i*(kx0-ky0)))*sigma_plus/5 ...
           + A(2,1)*sigma_z ...
           + A(1,1)*(1 + exp(-1i*(kx0+ky0)) + exp(-1i*(2*kx0)) + exp(-1i*(kx0-ky0)))*sigma_minus/5;
       
        % Higher-order Fourier components
        H1    = computeH(A(1,1), kx0, ky0,  1, sigma_plus, sigma_minus);
        Hfu1  = computeH(A(1,1), kx0, ky0, -1, sigma_plus, sigma_minus);
        H2    = computeH(A(1,1), kx0, ky0,  2, sigma_plus, sigma_minus);
        Hfu2  = computeH(A(1,1), kx0, ky0, -2, sigma_plus, sigma_minus);
        H3    = computeH(A(1,1), kx0, ky0,  3, sigma_plus, sigma_minus);
        Hfu3  = computeH(A(1,1), kx0, ky0, -3, sigma_plus, sigma_minus);
        H4    = computeH(A(1,1), kx0, ky0,  4, sigma_plus, sigma_minus);
        Hfu4  = computeH(A(1,1), kx0, ky0, -4, sigma_plus, sigma_minus);
        H5    = computeH(A(1,1), kx0, ky0,  5, sigma_plus, sigma_minus);
        Hfu5  = computeH(A(1,1), kx0, ky0, -5, sigma_plus, sigma_minus);
        
        % Build effective Hamiltonian
        Heff = buildHeff(H0, H1, Hfu1, H2, Hfu2, H3, Hfu3, H4, Hfu4, H5, Hfu5, omg, T);
        
        % Diagonalize
        [V, ~] = eigs(Heff, 2*(2*omg+1));
        worker_psi_1(ky, :) = V(:, 2*omg+2);
    end
    
    % Time evolution
    for ky = 1:N
        for t = 1:N
            sum1 = 0;
            sum2 = 0;
            for j = 1:2*omg+1
                phase = exp(1i*(j-omg-1)*freq*(t-1)/(N-1)*2*pi/freq);  % Modified time origin
                idx = 2*j-1;
                sum1 = sum1 + phase * worker_psi_1(ky, idx);
                sum2 = sum2 + phase * worker_psi_1(ky, idx+1);
            end
            worker_psi_1_t_up(ky, t) = sum1;
            worker_psi_1_t_down(ky, t) = sum2;
        end
    end
    
    % Store worker results
    psi_1_t_up(kx, :, :) = worker_psi_1_t_up;
    psi_1_t_down(kx, :, :) = worker_psi_1_t_down;
end

% Post-processing and visualization
sigma_x = [0, 1; 1, 0];
sigma_y = [0, -1i; 1i, 0];
sigma_z = [1, 0; 0, -1];

% Preallocate storage for points satisfying conditions
redPoints = [];
bluePoints = [];

theta = 0.0*pi;
phi = 4*pi/4;               % Equivalent to pi
phi1 = 3*pi/4;
theta1 = theta + pi;

t_start = floor(0.2*N);%nontrival 0.2N trival 1
t_end = floor(0.6*N);% nontrival 0.6N   trival N

for kx = 1:N
    for ky = 1:N
        for t = t_start:t_end
            psi1 = psi_1_t_up(kx, ky, t);
            psi2 = psi_1_t_down(kx, ky, t);
            delta = 0.02;   % Tolerance
            
            % Compute expectation values
            xx = real([conj(psi1), conj(psi2)] * sigma_x * [psi1; psi2]);
            yy = real([conj(psi1), conj(psi2)] * sigma_y * [psi1; psi2]);
            zz = real([conj(psi1), conj(psi2)] * sigma_z * [psi1; psi2]);
            
            % Check red point condition
            if abs(xx - sin(phi)*cos(theta)) < delta && ...
               abs(yy - sin(phi)*sin(theta)) < delta && ...
               abs(zz - cos(phi)) < delta
                redPoints = [redPoints; kx, ky, t];
            end
            
            % Check blue point condition
            if abs(xx - sin(phi1)*cos(theta1)) < delta && ...
               abs(yy - sin(phi1)*sin(theta1)) < delta && ...
               abs(zz - cos(phi1)) < delta
                bluePoints = [bluePoints; kx, ky, t];
            end
        end
    end
end

% Save data
%save('E:\mnist\trivallink.mat', 'redPoints', 'N', 'bluePoints');

% Plot results
figure;
if ~isempty(redPoints)
    plot3(redPoints(:,1)/N*2, redPoints(:,2)/N*2, redPoints(:,3)/N, 'r.', 'MarkerSize', 10);
    hold on;
end
if ~isempty(bluePoints)
    plot3(bluePoints(:,1)/N*2, bluePoints(:,2)/N*2, bluePoints(:,3)/N, 'b.', 'MarkerSize', 10);
end
xlabel('k_x');
ylabel('k_y');
zlabel('t');
grid on;
hold off;

% Optionally close parallel pool
% delete(gcp);

% ---------- Helper functions ----------

function H = computeH(A, kx0, ky0, n, sigma_plus, sigma_minus)
    phase_coeff = 1i/(2*pi*n) * A;
    term1 = (exp(-1i*2*pi/5*n) - 1) * (sigma_plus + sigma_minus);
    term2 = (exp(-1i*4*pi/5*n) - exp(-1i*2*pi/5*n)) * (exp(1i*(kx0-ky0))*sigma_plus + exp(-1i*(kx0-ky0))*sigma_minus);
    term3 = (exp(-1i*6*pi/5*n) - exp(-1i*4*pi/5*n)) * (exp(1i*(2*kx0))*sigma_plus + exp(-1i*(2*kx0))*sigma_minus);
    term4 = (exp(-1i*8*pi/5*n) - exp(-1i*6*pi/5*n)) * (exp(1i*(kx0+ky0))*sigma_plus + exp(-1i*(kx0+ky0))*sigma_minus);
    H = phase_coeff * (term1 + term2 + term3 + term4);
end

function Heff = buildHeff(H0, H1, Hfu1, H2, Hfu2, H3, Hfu3, H4, Hfu4, H5, Hfu5, omg, T)
    n_total = 2*(2*omg+1);
    Heff = zeros(n_total);
    base_energy = 2*omg;  % Energy offset
    
    % Diagonal blocks
    for j = 1:(2*omg+1)
        idx = (2*j-1):(2*j);
        Heff(idx, idx) = H0 + (j + base_energy) * (2*pi/T) * eye(2);
    end
    
    % Off-diagonal blocks (Fourier components)
    % n = ±1
    for j = 1:(2*omg)
        Heff(2*j-1, 2*j+1) = Hfu1(1,1);
        Heff(2*j,   2*j+2) = Hfu1(2,2);
        Heff(2*j-1, 2*j+2) = Hfu1(1,2);
        Heff(2*j,   2*j+1) = Hfu1(2,1);
        
        Heff(2*j+1, 2*j-1) = H1(1,1);
        Heff(2*j+2, 2*j)   = H1(2,2);
        Heff(2*j+2, 2*j-1) = H1(2,1);
        Heff(2*j+1, 2*j)   = H1(1,2);
    end
    
    % n = ±2
    for j = 1:(2*omg-1)
        Heff(2*j-1, 2*j+3) = Hfu2(1,1);
        Heff(2*j,   2*j+4) = Hfu2(2,2);
        Heff(2*j-1, 2*j+4) = Hfu2(1,2);
        Heff(2*j,   2*j+3) = Hfu2(2,1);
        
        Heff(2*j+3, 2*j-1) = H2(1,1);
        Heff(2*j+4, 2*j)   = H2(2,2);
        Heff(2*j+4, 2*j-1) = H2(2,1);
        Heff(2*j+3, 2*j)   = H2(1,2);
    end
    
    % n = ±3
    for j = 1:(2*omg-2)
        Heff(2*j-1, 2*j+5) = Hfu3(1,1);
        Heff(2*j,   2*j+6) = Hfu3(2,2);
        Heff(2*j-1, 2*j+6) = Hfu3(1,2);
        Heff(2*j,   2*j+5) = Hfu3(2,1);
        
        Heff(2*j+5, 2*j-1) = H3(1,1);
        Heff(2*j+6, 2*j)   = H3(2,2);
        Heff(2*j+6, 2*j-1) = H3(2,1);
        Heff(2*j+5, 2*j)   = H3(1,2);
    end
    
    % n = ±4
    for j = 1:(2*omg-3)
        Heff(2*j-1, 2*j+7) = Hfu4(1,1);
        Heff(2*j,   2*j+8) = Hfu4(2,2);
        Heff(2*j-1, 2*j+8) = Hfu4(1,2);
        Heff(2*j,   2*j+7) = Hfu4(2,1);
        
        Heff(2*j+7, 2*j-1) = H4(1,1);
        Heff(2*j+8, 2*j)   = H4(2,2);
        Heff(2*j+8, 2*j-1) = H4(2,1);
        Heff(2*j+7, 2*j)   = H4(1,2);
    end
    
    % n = ±5
    for j = 1:(2*omg-4)
        Heff(2*j-1, 2*j+9)  = Hfu5(1,1);
        Heff(2*j,   2*j+10) = Hfu5(2,2);
        Heff(2*j-1, 2*j+10) = Hfu5(1,2);
        Heff(2*j,   2*j+9)  = Hfu5(2,1);
        
        Heff(2*j+9, 2*j-1)  = H5(1,1);
        Heff(2*j+10, 2*j)   = H5(2,2);
        Heff(2*j+10, 2*j-1) = H5(2,1);
        Heff(2*j+9, 2*j)    = H5(1,2);
    end
end
