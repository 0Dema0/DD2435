%% Equilibrium constant
R = 8.314; % Gas constant, in Joule/(Kelvin*mol)
T = 293; % Tempertaure, in K
Delta_G = -3.5; % Standard Gibbs free energy, in kcal/mol
K = exp(-Delta_G * 4184 / (R * T)); % Equilibrium constant
fprintf('The equilibrium constant K is: %.4f\n', K);
% Comment: since K >> 1, it means that inside the cell we have a lot more
% substrate than product, and the cell is strongly drawn to the reaction (it looks for equilibrium)

%% Enzyme simulation
K_m = 8e-3; % Michaleis constant, in M (or mol/liter)
K_x = 3e-6; % Dissociation constant, in M (or mol/liter)
K_eq = 0.8e-3;
n = 2.5; % Hill coefficient
V_f = 1.8e-3; % Limiting rate, in M/s (mol/(liter*s))

S = 0:0.0001:0.05; % Concentration of F6P in M
X = 0:0.0000005:0.000005; % Concentration of FBP in M
figure;
for i = 1:length(X)
    j = Hill(K_m,K_x,K_eq,n,V_f,S,X(i)); % Calculate the reaction rate


    plot(S, j/V_f);
    xlabel('Concentration of F6P (M)');
    ylabel('Reaction Rate (j)');
    title('Reaction Rate vs Concentration of F6P for different concentration of FBP');
    xlabel('S [M]');
    ylabel('j/V_f');
    hold on;
end
legend(arrayfun(@(x) sprintf('X=%.1e', x), X, 'UniformOutput', false));
hold off;

%% Stochiometry
N = [1,-1,0;0,1,-1];

S0 = [1e-3;1e-6];
tspan = 0:0.1:2000;

options = odeset('RelTol',1e-8,'AbsTol',1e-12,'NonNegative',1:2);
j_in_values = [0.6e-6, 6e-6];
for k = 1:length(j_in_values)
    j_in = j_in_values(k);
    [t,y] = ode15s(@(t,S) fun(t,S,j_in), [0,2000], S0, options);
    figure;
    plot(t, y(:,1), '-o', t, y(:,2), '-x');
    legend('S1 (F6P)', 'S2 (FBP)', 'Location', 'best');
    xlabel('Time (s)');
    ylabel('Concentration (M)');
    title(sprintf('Dynamics with j_{in} = %.1e M/s', j_in));
end


%% Functions
function j = Hill(K_m,K_x,K_eq,n,V_f,S,X)
%HILL Computes reaction rate using Hill equation
%
%   j = Hill(S)
%
% Inputs:
%   S        - Concentration of F6P [M]
%   X        - Concentration of FBP [M]
%
% Output:
%   j        - Reaction rate [M/s]

    sigma_eq = K_eq/K_m;
    alpha = sigma_eq^(-n) - 1;

    sigma = S./K_m;
    csi = X/K_x;
    j = V_f * ((sigma.^n) ./ (sigma.^n + ((1+csi^n)/(1+alpha*csi^n))));
end

function dSdt = fun(t,S,j_in)
    N = [1,-1,0;0,1,-1];

    K_m = 8e-3; % Michaleis constant, in M
    K_x = 3e-6; % Dissociation constant, in M (or mol/liter)
    K_eq = 0.8e-3;
    n = 2.5; % Hill coefficient
    V_f = 1.8e-3; % Limiting rate, in M/s (mol/(liter*s))
    j_r = Hill(K_m,K_x,K_eq,n,V_f,S(1),S(2)); % Flux between, in M/s

    V = 60e-6; % Limiting rate, in M/s
    K_m = 10e-6; % Michaleis constant, in M
    j_out = V * (S(2)/(K_m+S(2))); % Flux out, in M/s

    j = [j_in;j_r;j_out];

    dSdt = N*j;
end