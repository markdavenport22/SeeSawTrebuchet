function [mcw, MaxStress, EnergyTotal] = AOptimise(l1,l2,h,mcw,s,m,w)
    format longG

    % Fixed values
    ar = 45; % Release Angle
    % Optmising values
    
    % Invented & technical variables
    g = 9.81;
    y = 90 - ar;
    as = 60; %FIXED
    P = l1.*sind(as);
    H = P + l1*sind(ar);
    E = 10.3*10^9; % Young's modulus
    p = 500; % Density
    % Minimum Velocity
    num1 = (g.*s^2);
    den1 = ((-2.*H)+2.*s.*tand(y)).*(cosd(y).^2);
    VelocityMin= sqrt(num1./den1);
    AngVelMin = VelocityMin ./ l1;
    %Rotational Energy
    I_m = m.*l1.^2;
    I_mcw = mcw.*l2.^2;
    I_beam = p.*h.*w.*(l1+l2).*(((h.^2 + w.^2)./(12))+((-l2+l1)./2).^2);
    SyRE = (1/2).*(I_m + I_mcw + I_beam).*(AngVelMin).^2;
    %Potential Energy
    CWPE = mcw.*g.*(l2*sind(as)+l2.*sind(ar));
    PrPE = m.*g.*(l1*sind(as)+l1.*sind(ar));
    L2PE = (h.*w.*(l2).*p).*g.*((1/2)*l2*sind(as)+(1/2)*l2*sind(ar));
    L1PE = (h.*w.*(l1).*p).*g.*((1/2)*l1*sind(as)+(1/2)*l1*sind(ar));
    %Stress
    Angles = -60:1:45;
    AngularVelocity = ((0:1:105)/105)*AngVelMin;
    AngularAcceleration = (AngVelMin.^2)./210;    
    L2_Fy = mcw.*AngularAcceleration.*l2 + mcw*g*cosd(Angles);
    Sig_L2_Surface = (6.*L2_Fy.*l2)/(w.*h.^2);    
    L2_Fx = mcw.*l2.*AngularVelocity.^2 + mcw.*g.*sind(Angles);
    Sig_L2_Axial = (L2_Fx)./(w.*h);    
    L1_Fy = m*g*cosd(Angles) - m.*AngularAcceleration.*l1;
    Sig_L1_Surface = abs((6.*L1_Fy.*l1)/(w.*h.^2));    
    L1_Fx =  m.*l1.*AngularVelocity.^2 - m.*g.*sind(Angles);
    Sig_L1_Axial = (L1_Fx)./(w.*h);    
    MaxStress = max([max(Sig_L2_Surface) max(Sig_L2_Axial) max(Sig_L1_Surface) max(Sig_L1_Axial)]);    
    % Beam Bending & Strain Energy
    num2 = 2.*(l1.^3).*(m*g*cosd(as)).^2;
    den2 = E.*w.*h.^3;
    SL1BB = num2/den2;    
    num3 = 2.*(l2.^3).*(mcw*g*cosd(as)).^2;
    den3 = E.*w.*h.^3;
    SL2BB = num3/den3;    
    num4 = 2.*(l2.^3).*(mcw*g*cosd(ar)).^2;
    den4 = E.*w.*h.^3;
    RL2BB = num4/den4;    
    num5 = (l1).*(m.*g.*sind(as)).^2;
    den5 = 2.*w.*h.*E;
    SL1SE = num5./den5;
    num6 = (l2).*(mcw.*g.*sind(as)).^2;
    den6 = 2.*w.*h.*E;
    SL2SE = num6./den6;
    num7 = (l2).*(mcw.*g.*sind(ar)).^2;
    den7 = 2.*w.*h.*E;
    RL2SE = num7./den7;
    
    %Equivilance
    EnergyIn = CWPE + L2PE + SL1BB + SL2BB + SL1SE + SL2SE;
    EnergyOut = SyRE + PrPE + L1PE + RL2BB + RL2SE;
    EnergyTotal = EnergyIn - EnergyOut;
end