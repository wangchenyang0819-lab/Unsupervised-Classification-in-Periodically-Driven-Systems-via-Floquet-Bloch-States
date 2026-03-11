clear

% Parameters
N = 41;                     % Number of discrete grid points
omg = 30;                   % Number of quasi-energy levels
T = 2*pi/(1*pi);            % Driving period
neng = 2*pi/T;              % Frequency
g = 1.6;                    % Coupling strength
A22 = 0.3*pi;               % Additional coupling parameter
A1 = 0.3*pi;                % Coupling coefficient 1
A2 = 0.1*pi;                % Coupling coefficient 2

% Initialize wave function arrays
bohanshu1 = zeros(N, 2*(2*omg+1));
bohanshut1 = zeros(N, N);
bohanshut2 = zeros(N, N);

% Pauli matrices
paolix = [0, 1; 1, 0];
paoliy = [0, -1i; 1i, 0];
paoliz = [1, 0; 0, -1];

% Main calculation loop
for k = 1:N
    k0 = (k-1)*2*pi/(N-1);
    
    % Define Hamiltonians
    H0 = A2*paoliy + (A1*paolix + A22*((cos(k0)*paolix) + (sin(k0))*paoliy));
    H1 = g*((cos(k0)*paolix) + (sin(k0))*paoliy);
    
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
        bohanshu1(k, j) = V(j, 2*omg+2);
    end
end

% Time evolution of wave functions
for k = 1:N
    for t = 1:N
        for j = 1:2*omg+1
            phase = exp(-1i*(j-omg-1)*(neng)*(t-1)/(N-1)*2*pi/neng);
            bohanshut1(k, t) = bohanshut1(k, t) + phase * bohanshu1(k, 2*j-1);
            bohanshut2(k, t) = bohanshut2(k, t) + phase * bohanshu1(k, 2*j);
        end
    end
end

% Calculate Bloch vector components
xx = zeros(N, N);
yy = zeros(N, N);
zz = zeros(N, N);

for k = 1:N
    for t = 1:N
        bohanshu = [bohanshut1(k, t); bohanshut2(k, t)];
        xx(k, t) = real(bohanshu' * paolix * bohanshu);
        yy(k, t) = real(bohanshu' * paoliy * bohanshu);
        zz(k, t) = real(bohanshu' * paoliz * bohanshu);
    end
end

% Save data
save('E:\mnist\kernelAIII31.mat', 'xx', 'yy', 'zz');