%% SYMBOLIC DIAGONALIZATION
% Symbolic environment MATLAB
syms U0 G dG
assume(U0, 'positive')

% Number of Brillouin zones
nBZ = 10;

% Reciprocal lattice vectors (symbolic)
G = sym(linspace(-nBZ*2*pi, nBZ*2*pi, 2*nBZ + 1));

% Matrix form of the effective Hamiltonian
N = length(G);
M = sym(zeros(N, N));
for i=1:N
    M(i,i) = 0.5*G(i)^2;
    if i < N
        M(i,i+1) = U0;
    end
    if i > 1
        M(i,i-1) = U0;
    end
end

disp(M);

% V eigenvectors matrix, D diagonal form associated with M
[V, D] = eig(M);
D = diag(D);            % Extracts the diagonal elements of D -> Eigenvalues




%% CASE U0 == 1
V1 = vpa(subs(V, U0, 1)); % Variable-Precision-Arithmetics evaluates symbolic
                          % expressions with a given precision d (here 32 digits)
E1 = vpa(subs(D, U0, 1)); 
                          
[E1, idx] = sort(E1);     % Sorted eigenvalues with relative permutation of indices (idx)
V1 = V1(:, idx);          % Sorted eigenvectors

GS_E1 = E1(1);            % Ground state
GS_eig1 = V1(:,1);        % Eigenvector ground state
fprintf('Ground state for U0 = 1:\n');
disp(GS_E1);
%disp(GS_eig1);

xx = linspace(-3, 3, 300);                  % Real-space grid
dx = (max(xx) - min(xx)) / (length(xx)-1);  % x Grid spacing
psi = zeros(size(xx));              

psi = zeros(size(xx));                      % Wavefunction
for k=1:length(GS_eig1)
    psi = psi + GS_eig1(k) * cos(G(k) * xx);
end
psi = psi/ sqrt(sum(abs(psi).^2) * dx);     % Normalization
psi = vpa(psi);   

plot(xx, abs(psi));  % Wave function
hold on
plot(xx, cos(2*pi*xx)); % Potential
grid on





%% CASE U0 == 100
V100 = vpa(subs(V, U0, 100)); % Variable-Precision-Arithmetics evaluates symbolic
                              % expressions with a given precision d (here 32 digits)
E100 = vpa(subs(D, U0, 100)); 
                          
[E100, idx] = sort(E100);     % Sorted eigenvalues with relative permutation of indices (idx)
V100 = V100(:, idx);          % Sorted eigenvectors

GS_E100 = E100(1);            % Ground state
GS_eig100 = V100(:,1);        % Eigenvector ground state
fprintf('Ground state for U0 = 100:\n');
disp(GS_E100);
%disp(GS_eig100);

xx = linspace(-3, 3, 400);                 % Real-space grid
dx = (max(xx) - min(xx)) / (length(xx)-1);  % x Grid spacing

psi = zeros(size(xx));                      % Wavefunction
for k=1:length(GS_eig100)
    psi = psi + GS_eig100(k) * cos(G(k) * xx);
end
psi = psi / sqrt(sum(abs(psi).^2) * dx);     % Normalization
psi = vpa(psi);   

plot(xx, abs(psi).^2);  % Wave function
hold on
%plot(xx, cos(2*pi*xx)); % Potential
grid on
