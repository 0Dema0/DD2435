% rng('shuffle');
% s = rng;
% save('rng_seed_new2.mat','s');

% load("rng_seed.mat",'s');
% rng(s);
load("rng_seed_new.mat",'s');
rng(s);

%% Values and tables
R = 8.314; % Gas constant, in Joule/(Kelvin*mol)
F = 96480; % Faradays constant, in Coulomb/mol
T = 293; % Temperature, in Kelvin

Ion = {'K+'; 'Na+'; 'Cl-'};
P     = [4.00; 0.12; 0.40]* 1e-9; %Permeability, in m/s
C_in  = [400; 50; 40]; % Intracellular concentration, in mM (millimolar)
C_out = [10; 460; 5]; % Extracellular concentration, in mM (millimolar)
z     = [1; 1; -1]; % Charge, dimensionless

IonTable = table(P, C_in, C_out, z, ...
    'RowNames', Ion);
disp(IonTable);

d_s = 100*1e-6; % Diameter of a neuron soma, in m
d_d = 1*1e-6; % Diameter of a dendritic spine head, in m
t_m = 100; % Membrane thickness, in Angstrom
c_m = 1e-2; % Membrane specific capacitance, in F/m^2
v_m = -50*1e-3; % Membrane potential, in V

%% 3: GHK equation
% Task 1:
V_rest = GHK_voltage(R, F, T, IonTable); % Potential, in V
disp(['Resting potential: ' num2str(V_rest*1e3, '%.2f') ' mV']);

% Task 2:
IonTable_temp = IonTable;   % Make a copy
[IonTable_temp.C_in, IonTable_temp.C_out] = deal(IonTable_temp.C_out, IonTable_temp.C_in);
V_rest_opposite = GHK_voltage(R, F, T, IonTable_temp); % Potential, in V
disp(['Resting potential (opposite): ' num2str(V_rest_opposite*1e3, '%.2f') ' mV']);
% We get the opposite sign

% Task 3:
IonTable_temp = IonTable;   % Make a copy
IonTable_temp{'Na+','P'} = 0;
IonTable_temp{'Cl-','P'} = 0;
V_rest_K = GHK_voltage(R, F, T, IonTable_temp); % Potential, in V
disp(['Potential with only potassium permeability nonzero: ' num2str(V_rest_K*1e3, '%.2f') ' mV']);
% We get a more negative number

% Task 4:
IonTable_temp = IonTable;   % Make a copy
IonTable_temp{'K+','P'} = 0;
IonTable_temp{'Cl-','P'} = 0;
V_rest_Na = GHK_voltage(R, F, T, IonTable_temp); % Potential, in V
disp(['Potential with only sodium permeability nonzero: ' num2str(V_rest_Na*1e3, '%.2f') ' mV']);
% We get a postive potential

% Task 5:
V1 = -70*1e-3; % Potential, in V
I1 = GHK_current(R, F, T, V1, IonTable); % Current density, in A/m^2
I1_table = table(I1*1e4, 'RowNames', IonTable.Properties.RowNames, ...
                 'VariableNames', {'Current density in A/cm^2'});
disp('Ionic current densities at V = -70 mV:');
disp(I1_table);
I1_tot = sum(I1); % Current density, in A/m^2
disp(['Total current density at V = -70 mV: ' num2str(I1_tot*1e4, '%.2f') ' A/cm^2']);

V2 = 0; % Potential, in V
I2 = GHK_current(R, F, T, V2, IonTable); % Current density, in A/m^2
I2_table = table(I2*1e4, 'RowNames', IonTable.Properties.RowNames, ...
                 'VariableNames', {'Current density in A/cm^2'});
disp('Ionic current densities at V = 0 mV:');
disp(I2_table);
I2_tot = sum(I2); % Current density, in A/m^2
disp(['Total current density at V = 0 mV: ' num2str(I2_tot*1e4, '%.2f') ' A/cm^2']);

% task 5; I(V) curve
V_vec = (-80:5:80) / 1000;

I_matrix = zeros(length(V_vec), height(IonTable));

for i = 1:length(V_vec)
    I_matrix(i,:) = GHK_current(R, F, T, V_vec(i), IonTable);  % Returns current density vectors for each ion
end

I_vec_tot = sum(I_matrix, 2);

figure;
plot(V_vec*1e3, I_matrix(:,1), 'r', 'DisplayName', 'K^+');
hold on;
plot(V_vec*1e3, I_matrix(:,2), 'b', 'DisplayName', 'Na^+');
plot(V_vec*1e3, I_matrix(:,3), 'g', 'DisplayName', 'Cl^-');
xlabel('Voltage [mV]'); ylabel('Current density [A/m^2]');
title('I–V curves');
legend; grid on;

figure;
semilogy(V_vec*1e3, abs(I_matrix(:,1)), 'r', 'DisplayName', 'K^+');
hold on;
semilogy(V_vec*1e3, abs(I_matrix(:,2)), 'b', 'DisplayName', 'Na^+');
semilogy(V_vec*1e3, abs(I_matrix(:,3)), 'g', 'DisplayName', 'Cl^-');
xlabel('Voltage [mV]'); ylabel('Current density [A/m^2] (log scale)');
title('I–V curves (log-scale)');
legend; grid on;

figure;
plot(V_vec*1e3, I_vec_tot, 'DisplayName', 'K^+ + Na^+ + Cl^-')
xlabel('Voltage [mV]'); ylabel('Current density [A/m^2]');
title("Total current density")
legend; grid on;


%% 4: A spherical membrane compartment
% task 6
I_d_s = -GHK_current(R, F, T, v_m, IonTable); % Current density, in A/m^2
I_d_s_table = table(I_d_s*1e4, 'RowNames', IonTable.Properties.RowNames, ...
    'VariableNames', {'Current density in A/cm^2'});
disp('Ionic current densities through neuron soma membrane at V = -50 mV:');
disp(I_d_s_table);
I_d_s_tot = sum(I_d_s); % Current density, in A/m^2
disp(['Total current density through neuron soma membrane at V = -50 mV: ' num2str(I_d_s_tot*1e4, '%.2f') ' A/cm^2']);

A_s = 4*pi*(d_s/2)^2; % Area of soma membrane, in m^2
I_s = I_d_s * A_s; % Current, in A
I_s_table = table(I_s, 'RowNames', IonTable.Properties.RowNames, ...
    'VariableNames', {'Current in A'});
disp('Ionic currents through neuron soma membrane at V = -50 mV:');
disp(I_s_table);
I_s_tot = sum(I_s); % Current, in A
disp(['Total current through neuron soma membrane at V = -50 mV: ' num2str(I_s_tot, '%.2e') ' A']);

% Task 7:
G_i = I_s./(v_m-V_rest); % Ion conductance, in A/V
I_G_i_table = table(G_i, 'RowNames', IonTable.Properties.RowNames, ...
    'VariableNames', {'Conductance in A/V'});
disp('Ionic conductances in neuron soma membrane at V = -50 mV:');
disp(I_G_i_table);
G_m = sum(G_i); % Ion conductance, in A/V
disp(['Total conductance of neuron soma membrane at V = -50 mV: ' num2str(G_m, '%.2e') ' A/V']);

% Task 8:
dt = 0.1e-3;
t_max = 50e-3;
time = 0:dt:t_max;
V_m = zeros(2, length(time));
I_inj = zeros(2, length(time));
V_m(:,1) = -50e-3;

for i = 1:length(time)-1
    for trace = 1:2
        I_s = -GHK_current(R, F, T, V_m(trace,i), IonTable);
        I_s_m = sum(I_s);

        if trace == 2
            if (i >= length(time)/5) && (i < length(time)/3)
                I_inj(trace,i) = 0.05;
            else
                I_inj(trace,i) = 0;
            end
        end

        % Euler forward
        V_m(trace,i+1) = V_m(trace,i) + (I_inj(trace,i) + I_s_m)/c_m * dt;
    end
end

for trace = 1:2
    figure;
    % Injected current
    subplot(2,1,1);
    plot(time*1e3, I_inj(trace,:)*1e3);
    xlabel('Time [ms]');
    ylabel('Injected current [mA]');
    ylim([-5,55]);
    title('I_{inj} vs time');

    subplot(2,1,2);
    plot(time*1e3, V_m(trace,:)*1e3);
    xlabel('Time [ms]');
    ylabel('Membrane potential [mV]');
    title('V_m vs time');

    if trace == 1
        sgtitle('Plots of I_{inj}(t) = 0 and V_m(t) vs Time (Soma)');
    else
        sgtitle('Plots of I_{inj}(t) ≠ 0 and V_m(t) vs Time (Soma)');
    end
end

% Task 9:
tau_m = zeros(1,2);
V_37 = zeros(1,2);
idx_max = zeros(1,2);
idx_37 = zeros(1,2);

for trace = 1:2
    V = V_m(trace,:);
    [maxV, idx] = max(V);
    idx_max(trace) = idx;
    V_37(trace) = (maxV - V_rest) * 0.37;
    [~, local_idx] = min(abs(V(idx:end) - (V_rest + V_37(trace))));
    idx_37(trace) = local_idx + idx - 1;
    tau_m(trace) = time(idx_37(trace)) - time(idx);
end

for trace = 1:2
    figure;
    
    subplot(2,1,1);
    plot(time*1e3, I_inj(trace,:)*1e3);
    xlabel('Time [ms]');
    ylabel('Injected current [mA]');
    ylim([-5,55]);
    title('I_{inj} vs time');
    
    subplot(2,1,2);
    plot(time*1e3, V_m(trace,:)*1e3);
    yline((V_rest + V_37(trace))*1e3, '--k', ...
        ['V_{37%} = ' num2str((V_rest + V_37(trace))*1e3,'%.2f') ' mV']);
    xline(time(idx_max(trace))*1e3, '--r', ...
        ['t_0 = ' num2str(time(idx_max(trace))*1e3,'%.2f') ' ms'], ...
        'LabelOrientation', 'horizontal');
    xline(time(idx_37(trace))*1e3, '--r', ...
        ['t_m = ' num2str(time(idx_37(trace))*1e3,'%.2f') ' ms'], ...
        'LabelOrientation', 'horizontal');
    xlabel('Time [ms]');
    ylabel('Membrane potential [mV]');
    title(['V_m vs time: \tau_m = t_m - t_0 = ' num2str(tau_m(trace)*1e3,'%.2f') 'ms']);
    
    if trace == 1
        sgtitle('Plots of I_{inj}(t) = 0 and V_m(t) vs Time with \tau_m and V_{37%} (Soma)');
    else
        sgtitle('Plots of I_{inj}(t) ≠ 0 and V_m(t) vs Time with \tau_m and V_{37%} (Soma)');
    end
end


% Task 10:
dt = 0.1e-3;
t_max = 50e-3;
time = 0:dt:t_max;
V_m = zeros(size(time));
V_m(1) = -50e-3;

for i = 1:length(time)-1
    IonTable_temp = IonTable;   % Make a copy
    if (0 <= i) && (i < 11) % From 0.000 to 0.010 seconds
        IonTable_temp{'K+','P'} = 4.00*1e-9;
        IonTable_temp{'Na+','P'} = 0.12*1e-9;
    elseif (11 <= i) && (i < 16) % From 0.010 to 0.015 seconds
        IonTable_temp{'K+','P'} = 4.00*1e-9;
        IonTable_temp{'Na+','P'} = 6.00*1e-9;
    elseif (16 <= i) && (i < 26) % From 0.015 to 0.025 seconds
        IonTable_temp{'K+','P'} = 4.00*1e-9;
        IonTable_temp{'Na+','P'} = 0.12*1e-9;
    elseif (26 <= i) && (i < 31) % From 0.025 to 0.030 seconds
        IonTable_temp{'K+','P'} = 40.0*1e-9;
        IonTable_temp{'Na+','P'} = 0.12*1e-9;
    else % From 0.030 onward
        IonTable_temp{'K+','P'} = 4.00*1e-9;
        IonTable_temp{'Na+','P'} = 0.12*1e-9;
    end
    I_s = -GHK_current(R,F,T,V_m(i),IonTable_temp);
    I_s_m = sum(I_s);

    I_inj = 0;

    % Euler forward
    V_m(i+1) = V_m(i) + (I_inj + I_s_m)/c_m * dt;
end

figure;
plot(time*1e3, V_m*1e3);
xlabel('Time [ms]');
ylabel('Membrane potential [mV]');
title('V_m course for soma with changing ion permeabilities (I_{inj} = 0)');

% Task 11:
A_d = 4*pi*(d_d/2)^2; % Area of dendritic spine head membrane, in m^2
I_d = I_d_s * A_d; % Current, in A
I_d_table = table(I_d, 'RowNames', IonTable.Properties.RowNames, ...
    'VariableNames', {'Current in A'});
disp('Ionic currents through neuron dendritic spine head membrane at V = -50 mV:');
disp(I_d_table);
I_d_tot = sum(I_d); % Current, in A
disp(['Total current through neuron dendritic spine head membrane at V = -50 mV: ' num2str(I_d_tot, '%.2e') ' A']);

%% 5: Stochastic, voltage dependent ion channels
% Task 12:
V_m = (-80:10:80)*1e-3;
m_inf = zeros(size(V_m));
tau_m_inf = zeros(size(V_m));

for i = 1:length(V_m)
    a_m = alpha_m(V_m(i));
    b_m = beta_m(V_m(i));
    m_inf(i) = a_m / (a_m + b_m);
    tau_m_inf(i) = 1 / (a_m + b_m);
end

figure;
plot(V_m*1e3, m_inf);
xlabel('V_m [mV]');
ylabel('m_{\infty} coefficient');
ylim([-0.1 1.1]);
title('m_{\infty} vs V');

figure;
plot(V_m*1e3, tau_m_inf);
xlabel('V_m [mV]');
ylabel('\tau_{m_{\infty}} coefficient');
title('\tau_{m_{\infty}} vs V');

% Task 13:
s_m = zeros(length(V_m),length(time)); % State of the Na-m channel (open or close)

for i = 1:length(V_m)
    a = alpha_m(V_m(i));
    b = beta_m(V_m(i));

    s_m(i,1) = rand() < m_inf(i);

    for j = 1:length(time)-1
        s_m(i,j+1) = next_state(s_m(i,j),a,b,dt);
    end
end

figure;
hold on;
for i = 1:length(V_m)
    plot(time*1e3, s_m(i,:));  % Time in ms, state 0 or 1
end
xlabel('Time [ms]');
ylabel('Na-m channel state (0=closed, 1=open)');
ylim([-0.1 1.1]);
title('Stochastic Na-m channel state vs time');
legend(arrayfun(@(V) sprintf('Vm = %.0f mV', V*1e3), V_m, 'UniformOutput', false));
hold off;

for i = 1:length(V_m)
    figure;
    plot(time*1e3, s_m(i,:));  % Time in ms, state 0 or 1
    xlabel('Time [ms]');
    ylabel(' Na-m channel state (0=closed, 1=open)');
    ylim([-0.1 1.1]);
    title(sprintf('Stochastic Na-m channel state at Vm = %.0f mV', V_m(i)*1e3));
end

% Task 14:
s_m1 = zeros(length(V_m), length(time)); % State of the first Na-m channel (open or close)
s_m2 = zeros(length(V_m), length(time)); % State of the second Na-m channel (open or close)
s_m3 = zeros(length(V_m), length(time)); % State of the third Na-m channel (open or close)
s_h = zeros(length(V_m),length(time)); % State of the Na-h channel (open or close)
s = zeros(length(V_m),length(time)); % State of the m^3h sodium channel (open or close)

for i = 1:length(V_m)
    a_m = alpha_m(V_m(i));
    b_m = beta_m(V_m(i));
    a_h = alpha_h(V_m(i));
    b_h = beta_h(V_m(i));

    s_m1(i,1) = rand() < m_inf(i);
    s_m2(i,1) = rand() < m_inf(i);
    s_m3(i,1) = rand() < m_inf(i);
    s_h(i,1) = rand() < a_h/(a_h + b_h);
    s(i,1) = (s_m(i,1)^3)*s_h(i,1);

    for j = 1:length(time)-1
        s_m1(i,j+1) = next_state(s_m1(i,j),a_m,b_m,dt);
        s_m2(i,j+1) = next_state(s_m2(i,j),a_m,b_m,dt);
        s_m3(i,j+1) = next_state(s_m3(i,j),a_m,b_m,dt);
        s_h(i,j+1) = next_state(s_h(i,j),a_h,b_h,dt);
        s(i,j+1) = s_m1(i,j+1)*s_m2(i,j+1)*s_m3(i,j+1)*s_h(i,j+1);
    end
end

figure;
hold on;
for i = 1:length(V_m)
    plot(time*1e3, s(i,:));  % Time in ms, state 0 or 1
end
xlabel('Time [ms]');
ylabel('Na (m^3h) channel state (0=closed, 1=open)');
ylim([-0.1 1.1]);
title('Stochastic Na (m^3h) channel state vs time');
legend(arrayfun(@(V) sprintf('Vm = %.0f mV', V*1e3), V_m, 'UniformOutput', false));
hold off;

for i = 1:length(V_m)
    figure;
    plot(time*1e3, s(i,:));  % Time in ms, state 0 or 1
    xlabel('Time [ms]');
    ylabel(' Na (m^3h) channel state (0=closed, 1=open)');
    ylim([-0.1 1.1]);
    title(sprintf('Stochastic Na (m^3h) channel state at Vm = %.0f mV', V_m(i)*1e3));
end

%% 6 Voluntary exercise
%Na_ch = 30;
Na_ch_vec = 0:10:30;
g_Na_max = 1e-12; % Max conductance for each Na channel, in S
E_Na = - (R*T/F) * log(IonTable{'Na+', 'C_out'}/IonTable{'Na+', 'C_in'});
C_m = c_m*A_d; % Example membrane capacitance in F
V_ch = -0.070 * ones(size(time)); % Initial Vm in V, e.g., -70 mV
V_leak = -0.070 * ones(size(time)); % Initial Vm in V, e.g., -70 mV
a_m_init = alpha_m(V_ch(1));
b_m_init = beta_m(V_ch(1));
a_h_init = alpha_h(V_ch(1));
b_h_init = beta_h(V_ch(1));
m_init = a_m_init / (a_m_init + b_m_init);
h_init = a_h_init / (a_h_init + b_h_init);

figure; hold on;

for Na_ch = Na_ch_vec
    I_Na = zeros(Na_ch,length(time));
    I_Na_tot = zeros(size(time));
    I_tot = zeros(size(time));
    I_leak = zeros(size(time));
    I_leak_inj = zeros(size(time));
    I_inj = zeros(size(time));
    
    s_m = zeros(3, Na_ch, length(time));
    s_h = zeros(Na_ch, length(time));
    s_ch = zeros(Na_ch, length(time));
    
    % Initialize channel states
    for ch = 1:Na_ch
        s_m(:,ch,1) = rand(3,1) < m_init;
        s_h(ch,1) = rand() < h_init;
    
        s_ch(ch,1) = prod(s_m(:,ch,1)) * s_h(ch,1);
        I_Na(ch,1) = g_Na_max * s_ch(ch,1) * (V_ch(1) - E_Na);
    end
    
    I_Na_tot(1) = sum(I_Na(ch,1));
    
    I_leak(1) = - sum(GHK_current(R,F,T,V_ch(1),IonTable)) * A_d;
    I_leak_inj(1) = I_leak(1);
    I_tot(1) = I_leak(1) + I_Na_tot(1);
    
    for j = 1:length(time)-1

        % Comment out next 4 lines for constant permeabilities
        IonTable_temp = IonTable;   % Make a copy
        if (Na_ch == 30) && (0.020 <= time(j+1) ) && (time(j+1) <= 0.025)
            IonTable_temp{'K+','P'} = 4.00*1e-8;
        end

        V_ch(j+1) = V_ch(j) + (I_tot(j)/C_m) * dt;
        V_leak(j+1) = V_leak(j) + (I_leak_inj(j)/C_m) * dt;
    
        for ch = 1:Na_ch
            for m = 1:3
                s_m(m,ch,j+1) = next_state(s_m(m,ch,j), alpha_m(V_ch(j)), beta_m(V_ch(j)), dt);
            end
            s_h(ch,j+1) = next_state(s_h(ch,j), alpha_h(V_ch(j)), beta_h(V_ch(j)), dt);
    
            s_ch(ch,j+1) = prod(s_m(:,ch,j+1)) * s_h(ch,j+1);
    
            I_Na(ch,j+1) = g_Na_max * s_ch(ch,j+1) * (V_ch(j) - E_Na);
        end
    
        if (0.010 <= time(j+1) ) && (time(j+1) <= 0.015)
            I_inj(j+1) = 4e-13;
            %I_inj(j+1) = 1e-13; Minimal no spike
        else
            I_inj(j+1) = 0;
        end
    
        I_Na_tot(j+1) = sum(I_Na(:,j+1),1);
    
        I_l1 = - GHK_current(R,F,T,V_leak(j+1),IonTable_temp);
        I_l2 = - GHK_current(R,F,T,V_ch(j+1),IonTable_temp);
    
        I_leak(j+1) = sum(I_l2) * A_d;
    
        I_leak_inj(j+1) = sum(I_l1) * A_d + I_inj(j+1);
        I_tot(j+1) = I_leak(j+1) + I_Na_tot(j+1) + I_inj(j+1);
    end
    
    plot(time*1e3, V_ch*1e3, 'LineWidth', 1.5, 'DisplayName', ...
        ['Na\_ch = ' num2str(Na_ch)]);

    if Na_ch == 30
        figure;
        subplot(2,1,1);
        plot(time*1e3, V_ch*1e3, 'LineWidth', 1.5);
        xlabel('Time [ms]');
        ylabel('V_{ch} [mV]');
        title('V_{ch} vs time');
        grid on;
        subplot(2,1,2);
        plot(time*1e3, V_leak*1e3, 'LineWidth', 1.5);
        xlabel('Time [ms]');
        ylabel('V_{leak} [mV]');
        title('V_{leak} vs time');
        grid on;
        sgtitle('Membrane potential with and without stochastic Na channels');

        figure;
        plot(time*1e3, I_Na_tot*1e12, 'LineWidth', 1.5);
        xlabel('Time [ms]');
        ylabel('Na current [pA]');
        title('Sodium current (I_{Na})');
        grid on;

        figure;
        plot(time*1e3, I_leak*1e12, 'LineWidth', 1.5);
        xlabel('Time [ms]');
        ylabel('Leak current [pA]');
        title('Leak current (I_{leak})');
        grid on;
    end
end

xlabel('Time [ms]');
ylabel('V_{ch} [mV]');
title('Membrane potential for different Na channel counts');
grid on;
hold off;

%% Functions
function V_rest = GHK_voltage(R, F, T, IonTable)
%GHK_VOLTAGE Compute resting potential using the GHK equation
%
%   V_rest = GHK_voltage(R, F, T, IonTable)
%
% Inputs:
%   R        - Gas constant [J/(K*mol)]
%   F        - Faraday constant [C/mol]
%   T        - Temperature [K]
%   IonTable - Table with row names {'K+','Na+','Cl-'} 
%              and variables: P, C_in, C_out
%
% Output:
%   V_rest   - Membrane potential [V]

    % Numerator (note Cl- uses C_in)
    num = IonTable{'K+','P'}  * IonTable{'K+','C_out'} + ...
          IonTable{'Na+','P'} * IonTable{'Na+','C_out'} + ...
          IonTable{'Cl-','P'} * IonTable{'Cl-','C_in'};

    % Denominator (note Cl- uses C_out)
    den = IonTable{'K+','P'}  * IonTable{'K+','C_in'} + ...
          IonTable{'Na+','P'} * IonTable{'Na+','C_in'} + ...
          IonTable{'Cl-','P'} * IonTable{'Cl-','C_out'};

    % GHK equation
    V_rest = (R*T/F) * log(num/den);
end

function I = GHK_current(R, F, T, V, IonTable)
%GHK_CURRENT Compute ionic current using GHK equation
%
%   I = GHK_current(R, F, T, V, IonTable)
%
% Inputs:
%   R        - Gas constant [J/(K*mol)]
%   F        - Faraday constant [C/mol]
%   T        - Temperature [K]
%   V        - Membrane potential [V]
%   IonTable - Table with row names {'K+','Na+','Cl-'}
%              and variables: P, C_in, C_out, z
%
% Output:
%   I        - Ionic current density [A/m^2] for each ion
%              (total current in Amperes = I * membrane area)

    % Extract table variables
    P     = IonTable.P;
    C_in  = IonTable.C_in;
    C_out = IonTable.C_out;
    z     = IonTable.z;

    % Dimensionless factor
    csi = z .* V * F / (R*T);

    % GHK equation
    I = zeros(height(IonTable));
    
    if V == 0
        I = P .* z .* F .* (C_in - C_out);
    else
        I = P .* z .* F .* csi .* ((C_in - C_out .* exp(-csi)) ./ (1 - exp(-csi)));
    end
end

function alpha_m = alpha_m(V)
%ALPHA Compute alpha_m value
%
%   alpha_m = alpha_m(V)
%
% Inputs:
%   V        - Membrane potential [V]
%
% Output:
%   alpha_m  - rate coefficient

    if V == -0.035
        alpha_m = 1e3;
    else
        alpha_m = -(V + 0.035)*1e5/(exp(-(V + 0.035)*1e2)-1);
    end
end

function beta_m = beta_m(V)
%ALPHA Compute beta_m value
%
%   beta_m = beta_m(V)
%
% Inputs:
%   V        - Membrane potential [V]
%
% Output:
%   beta_m   - rate coefficient

    beta_m = exp(-(V + 0.060)/0.018)*4e3;
end

function next_state = next_state(state,alpha,beta,dt)
    nch = size(state,2);
    p01 = rand(1,nch);
    alphadt = repmat(alpha,1,nch)*dt;
    betadt = repmat(beta,1,nch)*dt;
    next_state1 = (p01<alphadt) .* (state==0);
    next_state0 = (p01<betadt) .* (state==1);
    next_state = state + next_state1 - next_state0;
end

function alpha_h = alpha_h(V)
%ALPHA Compute alpha_h value
%
%   alpha_m = alpha_h(V)
%
% Inputs:
%   V        - Membrane potential [V]
%
% Output:
%   alpha_h  - rate coefficient

    alpha_h = 12*exp(-V/0.020);
end

function beta_h = beta_h(V)
%ALPHA Compute beta_h value
%
%   beta_m = beta_h(V)
%
% Inputs:
%   V        - Membrane potential [V]
%
% Output:
%   beta_h   - rate coefficient

    beta_h = 180/(exp(-(V + 0.030)*1e2) + 1);
end