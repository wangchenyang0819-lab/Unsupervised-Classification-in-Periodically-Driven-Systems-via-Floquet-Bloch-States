clear

% Pauli matrices
paolix = [0, 1; 1, 0];
paoliy = [0, -1i; 1i, 0];
paoliz = [1, 0; 0, -1];

% Physical parameters
omg = pi;                 % Driving frequency
N = 20;                   % Number of momentum points
T = 2*pi/omg;             % Driving period
NN = 600;                 % Parameter grid resolution
ww = 0.3*pi;              % Static coupling strength
g = 0.5*pi/T;             % Modulation strength

% Time discretization
time_steps = 60;          % Number of time steps
dt = T/time_steps;        % Time step size

% Initialize topological invariants arrays
raoshu00 = zeros(NN);
raoshupipi = zeros(NN);

% Calculate topological invariants for parameter grid
for a = 1:NN
    for b = 1:NN
        % Define parameters for current grid point
        v = 0.5*pi*(b/NN);
        w = 0.5*pi*((a-NN/2)/NN);
        
        % Initialize evolution operators
        UT = zeros(2, 2, N);       % Full period evolution
        UhalfT = zeros(2, 2, N);   % Half period evolution
        
        % Loop over momentum points
        for k_idx = 1:N
            k0 = 2*pi*(k_idx-1)/N;
            
            % Initialize evolution operators
            U_full = eye(2);
            U_half = eye(2);
            
            % Time evolution
            for t_step = 1:time_steps
                t_val = (t_step-1)*dt;
                
                % Time-dependent Hamiltonian
                H_t = (v + (ww + 2*g*cos(omg*t_val))*cos(k0))*paolix ...
                    + (ww + 2*g*cos(omg*t_val))*sin(k0)*paoliy ...
                    + w*paoliy;
                
                % Evolution step
                U_step = expm(-1i * H_t * dt);
                U_full = U_step * U_full;
                
                % Half period evolution
                if t_step <= time_steps/2
                    U_half = U_step * U_half;
                end
            end
            
            % Store evolution operators
            UT(:,:,k_idx) = U_full;
            UhalfT(:,:,k_idx) = U_half;
        end
        
        % Calculate effective Hamiltonians
        heff0 = zeros(2, 2, N);
        heffpi = zeros(2, 2, N);
        Ek = zeros(2, N);
        Vk = zeros(2, 2, N);
        angle0 = zeros(2, N);
        anglepi = zeros(2, N);
        Vrao0 = zeros(2, 2, N);
        Vraopi = zeros(2, 2, N);
        
        % Diagonalize full-period evolution operators
        for k = 1:N
            [V, E] = eig(UT(:,:,k));
            e = diag(E);
            Ek(:,k) = e;
            Vk(:,:,k) = V;
        end
        
        % Calculate quasi-energies
        for k = 1:N
            if imag(Ek(1,k)) <= 0
                anglepi(2,k) = -acos(real(Ek(1,k)))/T - 2*pi/T;
                anglepi(1,k) = acos(real(Ek(1,k)))/T - 2*pi/T;
            end
            if imag(Ek(1,k)) > 0
                anglepi(1,k) = -acos(real(Ek(1,k)))/T - 2*pi/T;
                anglepi(2,k) = acos(real(Ek(1,k)))/T - 2*pi/T;
            end
            
            if anglepi(1,k) <= -2*pi/T
                angle0(1,k) = anglepi(1,k) + 2*pi/T;
                angle0(2,k) = anglepi(2,k);
            end
            if anglepi(1,k) > -2*pi/T
                angle0(2,k) = anglepi(2,k) + 2*pi/T;
                angle0(1,k) = anglepi(1,k);
            end
        end
        
        % Construct effective Hamiltonians
        for k = 1:N
            heff0(:,:,k) = angle0(1,k) * (Vk(:,1,k) * Vk(:,1,k)') + ...
                          angle0(2,k) * (Vk(:,2,k) * Vk(:,2,k)');
            heffpi(:,:,k) = anglepi(1,k) * (Vk(:,1,k) * Vk(:,1,k)') + ...
                           anglepi(2,k) * (Vk(:,2,k) * Vk(:,2,k)');
        end
        
        % Construct winding operators
        for k = 1:N
            Vrao0(:,:,k) = UhalfT(:,:,k) * expm(1i*heff0(:,:,k)*T/2);
            Vraopi(:,:,k) = UhalfT(:,:,k) * expm(1i*heffpi(:,:,k)*T/2);
        end
        
        % Calculate winding numbers
        raoshu0 = 0;
        raoshupi = 0;
        for k = 1:N-1
            raoshu0 = raoshu0 + 1i/Vrao0(1,2,k)*(Vrao0(1,2,k+1)-Vrao0(1,2,k))/2/pi;
            raoshupi = raoshupi + 1i/Vraopi(1,1,k)*(Vraopi(1,1,k+1)-Vraopi(1,1,k))/2/pi;
        end
        
        % Round and store results
        raoshu00(a,b) = round(real(raoshu0));
        raoshupipi(a,b) = round(real(raoshupi));
        
        % Classify topological phases
        if raoshu00(a,b) == 1 && raoshupipi(a,b) == 0
            tpp(a,b) = 1;
        elseif raoshu00(a,b) == 1 && raoshupipi(a,b) == -1
            tpp(a,b) = 2;
        elseif raoshu00(a,b) == 0 && raoshupipi(a,b) == -1
            tpp(a,b) = 3;
        end
    end
end

% Plot phase diagram
figure;
axes1 = axes;
cmap = [0.5, 0.5, 0.5];
custom_colormap = [194/256, 206/256, 220/256; 
                   145/256, 173/256, 158/256; 
                   216/256, 156/256, 122/256; 
                   0.5, 0, 0.5];
x = 0:1/(NN-1):1;
y = -0.5:1/(NN-1):0.5;
colormap(custom_colormap);
caxis([1, 4]);
hold on;
surf(x, y, tpp, 'FaceColor', cmap);
ax = gca;
shading interp;
set(ax, 'TickLabelInterpreter', 'latex', ...
    'XTick', [0, 0.5, 1], 'XTickLabel', {'0', '0.5', '1'}, ...
    'YTick', [-0.5, 0, 0.5], 'YTickLabel', {'-0.5', '0', '0.5'}, ...
    'FontSize', 32);
xlabel('$\gamma_2 T/\pi$', 'Interpreter', 'latex', 'Color', 'k');
ylabel('$\gamma_3 T/\pi$', 'Interpreter', 'latex', 'Color', 'k');

% Second part: Cluster analysis
fangcha = 0.5;
l = 10;
hen = 16;
m = hen*hen;                    % Number of samples
A = zeros(2, m);                % Coupling coefficients
N = 60;                         % Number of momentum points
omg = 10;                       % Number of quasi-energy levels
T = 2*pi/pi;                    % Driving period
neng = 2*pi/T;                  % Frequency
g = 1.6;                        % Coupling strength
A22 = 0.3*pi;                   % Additional coupling

% Set parameter grid
for i = 1:hen
    for j = 1:hen
        A(1, (i-1)*hen+j) = 0.5*pi*((j-0.5)/(hen+0.5));
        A(2, (i-1)*hen+j) = 0.5*pi*((i-hen/2-0.5)/(hen+0.5));
    end
end

% Initialize arrays
bohanshu1 = zeros(N, 2*(2*omg+1), m);
bohanshu2 = zeros(N, 2*(2*omg+1), m);
bohanshut1 = zeros(N, N, m);
bohanshut2 = zeros(N, N, m);
bohanshut11 = zeros(N, N, m);
bohanshut22 = zeros(N, N, m);
Heff = zeros(2*(2*omg+1), 2*(2*omg+1));
z = zeros(1, m);
neihe = zeros(m, m);
touying0 = zeros(2, 2, N, N, m);

% Main calculation loop
for i = 1:m
    for k = 1:N
        k0 = k*2*pi/N;
        
        % Define Hamiltonians
        H0 = A(2,i)*paoliy + (A(1,i)*paolix + A22*((cos(k0)*paolix) + (sin(k0))*paoliy));
        H1 = g*((cos(k0)*paolix) + (sin(k0))*paoliy);
        Hfu1 = H1;
        
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
        
        % Diagonalize
        [V, E] = eigs(Heff, 2*(2*omg+1));
        e1 = diag(real(E) - (2*omg+(omg+1)*2*pi/T)*eye(2*(2*omg+1)));
        
        for j = 1:2*(2*omg+1)
            bohanshu1(k,j,i) = V(j, 2*omg+2);
            bohanshu2(k,j,i) = V(j, 2*omg+1);
        end
    end
end

% Time evolution
for i = 1:m
    for k = 1:N
        for t = 1:N
            for j = 1:2*omg+1
                phase = exp(-1i*(j-omg-1)*(neng)*t/N*2*pi/neng);
                bohanshut1(k,t,i) = bohanshut1(k,t,i) + phase * bohanshu1(k,2*j-1,i);
                bohanshut2(k,t,i) = bohanshut2(k,t,i) + phase * bohanshu1(k,2*j,i);
                bohanshut11(k,t,i) = bohanshut11(k,t,i) + phase * bohanshu2(k,2*j-1,i);
                bohanshut22(k,t,i) = bohanshut22(k,t,i) + phase * bohanshu2(k,2*j,i);
            end
        end
    end
end

% Calculate projection operators
for i = 1:m
    for kx = 1:N
        for t = 1:N
            psi1 = [bohanshut1(kx,t,i); bohanshut2(kx,t,i)];
            psi2 = [bohanshut11(kx,t,i); bohanshut22(kx,t,i)];
            touying0(:,:,kx,t,i) = 1 * (psi1 * psi1') - 1 * (psi2 * psi2');
        end
    end
end

% Calculate similarity matrix
for j = 1:m
    for o = 1:m
        neihe(o,j) = 1;
        for kx = 1:N
            for t = 1:N
                det_val = det(touying0(:,:,kx,t,o) + touying0(:,:,kx,t,j));
                neihe(o,j) = (1 - exp(-abs(det_val)^2/(fangcha^2))) * neihe(o,j);
            end
        end
    end
end

% Normalize and perform clustering
for i = 1:m
    z(1,i) = sum(neihe(i,:));
end

P = zeros(m, m);
for i = 1:m
    for j = 1:m
        P(i,j) = neihe(i,j) / sqrt(z(1,i)) / sqrt(z(1,j));
    end
end

% Eigenvalue decomposition
[V, E] = eigs(P, m);
e = diag(E);

% Hierarchical clustering
dataa = V(:, 1:3);
data = dataa;
ZZ = linkage(data, 'single');
TT = cluster(ZZ, 'maxclust', 3);

% Create phase map
f = zeros(hen, hen);
for i = 1:hen
    for j = 1:hen
        for m_cluster = 1:25
            if TT((i-1)*hen+j,1) == m_cluster
                f(i,j) = m_cluster;
            end
        end
    end
end

% Create parameter arrays
for i = 1:hen*hen
    B(1,i) = A(1,i)/pi*T;
    B(2,i) = A(2,i)/pi*T;
end

% Plot clusters
for i = 1:hen
    for j = 1:hen
        if f(i,j) == 1
            plot3(B(1,(i-1)*hen+j), B(2,(i-1)*hen+j), 10, '.', 'color', [1 0 0], 'MarkerSize', 40);
            hold on;
        elseif f(i,j) == 2
            plot3(B(1,(i-1)*hen+j), B(2,(i-1)*hen+j), 10, '.', 'color', [0 1 0], 'MarkerSize', 40);
            hold on;
        elseif f(i,j) == 3
            plot3(B(1,(i-1)*hen+j), B(2,(i-1)*hen+j), 10, '.', 'color', [0 0 1], 'MarkerSize', 40);
            hold on;
        end
    end
end

hold on;
set(gca, 'FontName', 'Times New Roman', 'fontsize', 32);
axis square;
set(gca, 'LineWidth', 5);
box on;
set(gca, 'XColor', [0.5, 0.5, 0.5], 'YColor', [0.5, 0.5, 0.5], 'Layer', 'top');
ax = gca;
set(ax, 'TickLabelInterpreter', 'latex', ...
    'XTick', [0, 0.5, 1], 'XTickLabel', {'0', '0.5', '1'}, ...
    'YTick', [-0.5, 0, 0.5], 'YTickLabel', {'-0.5', '0', '0.5'}, ...
    'FontSize', 32);
ax.XAxis.TickLabelColor = 'k';
ax.YAxis.TickLabelColor = 'k';
xlabel('$\gamma_2 T/\pi$', 'Interpreter', 'latex', 'Color', 'k');
ylabel('$\gamma_3 T/\pi$', 'Interpreter', 'latex', 'Color', 'k');

% PCA analysis
X = V(:, 1:3);
[coeff, score, latent, tsquared, explained, mu] = pca(X, 'Algorithm', 'svd');
figure;
plot(score(:,1), score(:,2), 'r*', 'MarkerSize', 10);

% Save data
save('E:\mnist\AIII.mat', 'e', 'score', 'tpp', 'B', 'f');