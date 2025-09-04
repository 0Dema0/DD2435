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

d_s = 100*1e-9; % Diameter of a neuron soma, in m
d_d = 1*1e-9; % Diameter of a dendritic spine head, in m
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
V_all = {V_m_1, V_m_2};
I_all = {I_inj_1, I_inj_2};
tau_m = zeros(1,2);
V_37 = zeros(1,2);
idx_max = zeros(1,2);
idx_37 = zeros(1,2);

for trace = 1:2
    V = V_all{trace};
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
    plot(time*1e3, I_all{trace}*1e3);
    xlabel('Time [ms]');
    ylabel('Injected current [mA]');
    ylim([-5,55]);
    title('I_{inj} vs time');
    
    subplot(2,1,2);
    plot(time*1e3, V_all{trace}*1e3);
    yline((V_rest + V_37(trace))*1e3, '--g', ...
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

%Task 11:
A_d = 4*pi*(d_d/2)^2; % Area of soma membrane, in m^2
I_d = I_d_s * A_d; % Current, in A
I_d_table = table(I_d, 'RowNames', IonTable.Properties.RowNames, ...
    'VariableNames', {'Current in A'});
disp('Ionic currents through neuron dendritic spine head membrane at V = -50 mV:');
disp(I_d_table);
I_d_tot = sum(I_d); % Current, in A
disp(['Total current through neuron dendritic spine head membrane at V = -50 mV: ' num2str(I_d_tot, '%.2e') ' A']);

%% 5: Stochastic, voltage dependent ion channels
%Task 12:
