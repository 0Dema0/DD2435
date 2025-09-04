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
