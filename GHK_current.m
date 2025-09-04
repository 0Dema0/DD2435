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


