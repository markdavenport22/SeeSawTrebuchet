function [Mass_Total, MaxAP, MaxBP] = BOptimise(OAP, OBP, d, t)
    format longG
    % Value Treatment Section
    OBP = 180 - OBP;
    disp("OAP:" + OAP + " OBP:" + OBP + " d:" + d + " t:" + t)
    
    % Values from Task A Decleration Section
    m = 15;
    l1 = 12.20; % L1 Length
    l2 = 2.35; % L2 Length
    mcw = 1037; % Counterweight Mass
    AngVelMax = 2.44365315; % The max value in the operation
    P = 10.565509926; % Insert P from optimisation task

    % Technical Values Decleration Section
    E = 200*10^9; % The youngs modlus of steel
    p = 7850; % density of steel

    % Stiffness Matrix Section
    % Generates stiffness matrices for the individual elements and the
    % system overall
    OAP_A = cosd(OAP); % A represents lambda
    OAP_U = sind(OAP); % U represents mu
    OBP_A = cosd(OBP); % A represents lambda
    OBP_U = sind(OBP); % U represents mu
    LAP = (P)/(sind(OAP));
    LBP = (P)/(sind(OBP));
    A = pi*(d/2)^2 - pi*(d/2 - t)^2;
    KAP = ((A*E)/LAP)*[(OAP_A^2) (OAP_A*OAP_U) 0 0 (-OAP_A^2) (-OAP_A*OAP_U);(OAP_A*OAP_U) (OAP_U^2) 0 0 (-OAP_A*OAP_U) (-OAP_U^2); 0 0 0 0 0 0 ; 0 0 0 0 0 0 ;(-OAP_A^2) (-OAP_A*OAP_U) 0 0 (OAP_A^2) (OAP_A*OAP_U); (-OAP_A*OAP_U) (-OAP_U^2) 0 0 (OAP_A*OAP_U) (OAP_U^2)];
    KBP = ((A*E)/LBP)*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 (OBP_A^2) (OBP_A*OBP_U) (-OBP_A^2) (-OBP_A*OBP_U); 0 0 (OBP_A*OBP_U) (OBP_U^2) (-OBP_A*OBP_U) (-OBP_U^2); 0 0 (-OBP_A^2) (-OBP_A*OBP_U) (OBP_A^2) (OBP_A*OBP_U); 0 0 (-OBP_A*OBP_U) (-OBP_U^2) (OBP_A*OBP_U) (OBP_U^2)];
    KGlobal = KAP + KBP;
    
    % Forces & Velocities Section
    % The forces and associated velocities are determined w.r.t. the angle
    % of the beam. 
    Angle = -60:1:45;
    w = ((Angle + 60)./105).*AngVelMax;
    F = mcw.*l2.*w.^2 - m*l1*w.^2;
    Fy = F.*sind(Angle) - mcw * 9.81 - m*9.81;
    Fx = F.*cosd(Angle);
    UC = (1/KGlobal(5,5)).*Fx;
    VC  = (1/KGlobal(5,5)).*Fy;
    
    % Stresses Section
    % Generates the stresses over time of each of the elements w.r.t the
    % angle of the beam. Also determines the max value achieved in either
    % beams to ensure dosen't reach yield strength.
    KAPe = [KAP(5,5) KAP(5,6) ; KAP(5,6) KAP(6,6)];
    KBPe = [KBP(5,5) KBP(5,6) ; KBP(5,6) KBP(6,6)];
    Sig_AP = (1/A)*[OAP_A OAP_U]*KAPe*[UC;VC];
    MaxAP = max(abs(Sig_AP));
    Sig_BP = (1/A)*[OBP_A OBP_U]*KBPe*[UC;VC];
    MaxBP = max(abs(Sig_BP));
    
    Mass_Total = p*A*(LAP + LBP);
end