clear

% Parameters
l = 10;                     % Dimension for phase diagram
m = 1;                      % Number of samples
A = zeros(2, m);            % Coupling coefficients
N = 200;                    % Number of discrete grid points
omg = 10;                   % Number of quasi-energy levels
T = 2*pi;                   % Driving period
neng = 2*pi/T;              % Frequency
g = 1*pi/T;                 % Coupling parameter

% Initialize parallel pool
if isempty(gcp('nocreate'))
    parpool;  % Start parallel pool with default settings
end

% Preallocate large arrays
bohanshu1 = zeros(N, N, 2*(2*omg+1));
bohanshut1 = zeros(N, N, N);
bohanshut2 = zeros(N, N, N);

% Set sample parameters
f = 4;%27
A(1,1) = (3*((f)/30))*pi/T;  % Hopping amplitude
A(2,1) = 0.5*pi/T;           % Sublattice potential

% Parallel loop over momentum points
parfor kx = 1:N
    % Each worker has its own variable copies
    paolix = [0,1;1,0];
    paoliy = [0,-1i;1i,0];
    paoliz = [1,0;0,-1];
    
    % Preallocate worker-local variables
    worker_bohanshu1 = zeros(N, 2*(2*omg+1));
    worker_bohanshut1 = zeros(N, N);
    worker_bohanshut2 = zeros(N, N);
    
    kx0 = sqrt(2)*(kx-1-N/2)*2*pi/(N-1);
    
    for ky = 1:N
        ky0 = sqrt(2)*(ky-1-N/2)*2*pi/(N-1);
        
        % Define static Hamiltonian
        H0 = 2/5*A(1,1)*(cos(kx0)+cos(ky0))*paolix + A(2,1)*paoliz;
        paolizheng = (paolix+1i*paoliy)/2;
        paolifu = (paolix-1i*paoliy)/2;
        
        % Compute Fourier components of driving Hamiltonian
        H1 = computeH(A(1,1), kx0, ky0, 1, paolizheng, paolifu);
        Hfu1 = computeH(A(1,1), kx0, ky0, -1, paolizheng, paolifu);
        H2 = computeH(A(1,1), kx0, ky0, 2, paolizheng, paolifu);
        Hfu2 = computeH(A(1,1), kx0, ky0, -2, paolizheng, paolifu);
        H3 = computeH(A(1,1), kx0, ky0, 3, paolizheng, paolifu);
        Hfu3 = computeH(A(1,1), kx0, ky0, -3, paolizheng, paolifu);
        H4 = computeH(A(1,1), kx0, ky0, 4, paolizheng, paolifu);
        Hfu4 = computeH(A(1,1), kx0, ky0, -4, paolizheng, paolifu);
        H5 = computeH(A(1,1), kx0, ky0, 5, paolizheng, paolifu);
        Hfu5 = computeH(A(1,1), kx0, ky0, -5, paolizheng, paolifu);
        
        % Build effective (Floquet) Hamiltonian
        Heff = buildHeff(H0, H1, Hfu1, H2, Hfu2, H3, Hfu3, H4, Hfu4, H5, Hfu5, omg, T);
        
        % Diagonalize effective Hamiltonian
        [V, ~] = eigs(Heff, 2*(2*omg+1));
        worker_bohanshu1(ky, :) = V(:, 2*omg+2);
    end
    
    % Time evolution (Fourier synthesis)
    for ky = 1:N
        for t = 1:N
            sum1 = 0;
            sum2 = 0;
            for j = 1:2*omg+1
                phase = exp(1i*(j-omg-1)*neng*(t-1)/(N-1)*2*pi/neng);
                idx = 2*j-1;
                sum1 = sum1 + phase * worker_bohanshu1(ky, idx);
                sum2 = sum2 + phase * worker_bohanshu1(ky, idx+1);
            end
            worker_bohanshut1(ky, t) = sum1;
            worker_bohanshut2(ky, t) = sum2;
        end
    end
    
    % Store results from worker
    bohanshut1(kx, :, :) = worker_bohanshut1;
    bohanshut2(kx, :, :) = worker_bohanshut2;
end

% Post-processing: calculate expectation values
paolix = [0,1;1,0];
paoliy = [0,-1i;1i,0];
paoliz = [1,0;0,-1];

% Arrays to store points that satisfy certain conditions
redPoints = [];
bluePoints = [];

% Target values for Bloch vector components
theta = 0.*pi;
phi = 3*pi/4;
phi1 = 3*pi/4;
theta1 = theta+pi;

% Time window for analysis
t_start = floor(1);
t_end = floor(1*N);

% Search for points that match target spin configurations
for kx = 1:N
    for ky = 1:N
        for t = t_start:t_end
            psi1 = bohanshut1(kx, ky, t);
            psi2 = bohanshut2(kx, ky, t);
            delta = 0.03;  % Tolerance for matching
            
            % Calculate Bloch vector components
            xx = real([conj(psi1), conj(psi2)] * paolix * [psi1; psi2]);
            yy = real([conj(psi1), conj(psi2)] * paoliy * [psi1; psi2]);
            zz = real([conj(psi1), conj(psi2)] * paoliz * [psi1; psi2]);
            
            % Check red points condition
            if abs(xx - sin(phi)*cos(theta)) < delta && ...
               abs(yy - sin(phi)*sin(theta)) < delta && ...
               abs(zz - cos(phi)) < delta
                redPoints = [redPoints; kx, ky, t];
            end
            
            % Check blue points condition
            if abs(xx - sin(phi1)*cos(theta1)) < delta && ...
               abs(yy - sin(phi1)*sin(theta1)) < delta && ...
               abs(zz - cos(phi1)) < delta
                bluePoints = [bluePoints; kx, ky, t];
            end
        end
    end
end

% Visualize the results
figure;
if ~isempty(redPoints)
    plot3(redPoints(:,1)/N*2, redPoints(:,2)/N*2, redPoints(:,3)/N, 'r.', 'MarkerSize', 15);
    hold on;
end
if ~isempty(bluePoints)
    plot3(bluePoints(:,1)/N*2, bluePoints(:,2)/N*2, bluePoints(:,3)/N, 'b.', 'MarkerSize', 15);
end
xlabel('kx');
ylabel('ky');
zlabel('t');
grid on;
hold off;

% Helper function: compute Fourier component of driving Hamiltonian
function H = computeH(A, kx0, ky0, n, paolizheng, paolifu)
    phase_coeff = 1i/(2*pi*n) * A;
    term1 = (exp(-1i*2*pi/5*n)-1) * (exp(1i*kx0)*paolizheng + exp(-1i*kx0)*paolifu);
    term2 = (exp(-1i*4*pi/5*n)-exp(-1i*2*pi/5*n)) * (exp(1i*ky0)*paolizheng + exp(-1i*ky0)*paolifu);
    term3 = (exp(-1i*6*pi/5*n)-exp(-1i*4*pi/5*n)) * (exp(-1i*kx0)*paolizheng + exp(1i*kx0)*paolifu);
    term4 = (exp(-1i*8*pi/5*n)-exp(-1i*6*pi/5*n)) * (exp(-1i*ky0)*paolizheng + exp(1i*ky0)*paolifu);
    H = phase_coeff * (term1 + term2 + term3 + term4);
end

% Helper function: build effective Hamiltonian in extended Hilbert space
function Heff = buildHeff(H0, H1, Hfu1, H2, Hfu2, H3, Hfu3, H4, Hfu4, H5, Hfu5, omg, T)
    n = 2*(2*omg+1);
    Heff = zeros(n);
    base_energy = 2*omg;  % Energy offset
    
    % Diagonal blocks
    for j = 1:(2*omg+1)
        idx = 2*j-1:2*j;
        Heff(idx, idx) = H0 + (j + base_energy) * (2*pi/T) * eye(2);
    end
    
    % Off-diagonal blocks (n=±1 Fourier components)
    for j = 1:(2*omg)
        Heff(2*j-1, 2*j+1) = Hfu1(1,1);
        Heff(2*j, 2*j+2) = Hfu1(2,2);
        Heff(2*j-1, 2*j+2) = Hfu1(1,2);
        Heff(2*j, 2*j+1) = Hfu1(2,1);
        
        Heff(2*j+1, 2*j-1) = H1(1,1);
        Heff(2*j+2, 2*j) = H1(2,2);
        Heff(2*j+2, 2*j-1) = H1(2,1);
        Heff(2*j+1, 2*j) = H1(1,2);
    end
    
    % n=±2 Fourier components
    for j = 1:(2*omg-1)
        Heff(2*j-1, 2*j+3) = Hfu2(1,1);
        Heff(2*j, 2*j+4) = Hfu2(2,2);
        Heff(2*j-1, 2*j+4) = Hfu2(1,2);
        Heff(2*j, 2*j+3) = Hfu2(2,1);
        
        Heff(2*j+3, 2*j-1) = H2(1,1);
        Heff(2*j+4, 2*j) = H2(2,2);
        Heff(2*j+4, 2*j-1) = H2(2,1);
        Heff(2*j+3, 2*j) = H2(1,2);
    end
    
    % n=±3 Fourier components
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
    
    % n=±4 Fourier components
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
    
    % n=±5 Fourier components
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
end
 % save('E:\mnist\trivallink.mat','redPoints','N','bluePoints');%xxzclassicalsmall
