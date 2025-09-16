%% Equilibrium constant
R = 8.314; % Gas constant, in Joule/(Kelvin*mol)
T = 293; % Tempertaure, in K
Delta_G = -3.5; % Standard Gibbs free energy, in kcal/mol
K = exp(-Delta_G * 4184 / (R * T)); % Equilibrium constant
fprintf('The equilibrium constant K is: %.4f\n', K);
% Comment: since K >> 1, it means that inside the cell we have a lot more
% substrate than prduct, and the cell is strongly drawn to the reaction (it look for equilibrium)

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
    ylabel('j/V_f []');
    hold on;
end
legend;
hold off;

%% Stochiometry
N = [1,-1,0;0,1,-1];

S0 = [1e-3;1e-9];

[t,y] = ode15s(@fun,[0,2000],S0);
figure;
plot(t,y(:,2),'-o');


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
%   j        - Reaction rate []

    sigma_eq = K_eq/K_m;
    alpha = sigma_eq^(-n) - 1;

    sigma = S./K_m;
    csi = X/K_x;
    j = V_f * ((sigma.^n) ./ (sigma.^n + ((1+csi^n)/(1+alpha*csi^n))));
end

function dSdt = fun(t,S)
    N = [1,-1,0;0,1,-1];

    j_in = 6e-6;

    K_m = 8e-3;
    K_x = 3e-6; % Dissociation constant, in M (or mol/liter)
    K_eq = 0.8e-3;
    n = 2.5; % Hill coefficient
    V_f = 1.8e-3; % Limiting rate, in M/s (mol/(liter*s))
    j_r = Hill(K_m,K_x,K_eq,n,V_f,S(1),S(2));

    V = 60e-6;
    K_m = 10e-6;
    j_out = V * (S(2)/(K_m+S(2)));

    j = [j_in;j_r;j_out];

    dSdt = N*j;
end


% T = 0:1:2000;
% S = zeros(size(T));
% X = zeros(size(T));
% j_in = zeros(size(T));
% j = zeros(size(T));
% j_out = zeros(size(T));
% 
% S(1) = 1e-3;
% X(1) = 1e-6;
% j_in(1) = 0.6e-6;
% j(1) = Hill(S(1),X(1));
% V = 60e-6;
% K_m = 10e-6;
% j_out(1) = V * (X(1)/(K_m+X(1)));