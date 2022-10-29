function [ OopDefl1, IPDefl1, TipDzb1, TTDspFA, TTDspSS, LSSTipV, NacYaw, TwrClrnB1, TwrClrnB2, TwrClrnB3, BldTwrHrClr, Surge, Sway, Heave, Roll, Pitch, Yaw] = Results( t,q, Controls, DOFs, ElastoDyn, Twr, Bld )

TipRad = ElastoDyn.TipRad;
HubRad = ElastoDyn.HubRad;
PreCone(1) = ElastoDyn.PreCone(1);
PreCone(2) = ElastoDyn.PreCone(2);
PreCone(3) = ElastoDyn.PreCone(3);
HubCM = ElastoDyn.HubCM;
OverHang = ElastoDyn.OverHang;
ShftTilt = ElastoDyn.ShftTilt;
NacCMxn = ElastoDyn.NacCMxn;
NacCMyn = ElastoDyn.NacCMyn;
NacCMzn = ElastoDyn.NacCMzn;
Twr2Shft = ElastoDyn.Twr2Shft;
TwrHt = ElastoDyn.TwrHt;
TowerBsHt = ElastoDyn.TowerBsHt;  
PtfmCMzt = ElastoDyn.PtfmCMzt;   
PtfmRefzt = ElastoDyn.PtfmRefzt;   

nt = Twr.nt;
TwrSec = Twr.TwrSec;
O1_TFA = Twr.O1_TFA;
O2_TFA = Twr.O2_TFA; 
O1_TSS = Twr.O1_TSS;
O2_TSS = Twr.O2_TSS;
s11_TFA = Twr.s11_TFA;
s12_TFA = Twr.s12_TFA;
s22_TFA = Twr.s22_TFA;
S11_TFA = Twr.S11_TFA;
S12_TFA = Twr.S12_TFA;
S22_TFA = Twr.S22_TFA;
s11_TSS = Twr.s11_TSS;
s12_TSS = Twr.s12_TSS;
s22_TSS = Twr.s22_TSS;
S11_TSS = Twr.S11_TSS;
S12_TSS = Twr.S12_TSS;
S22_TSS = Twr.S22_TSS;
dO1_TFA = Twr.dO1_TFA;
dO1_TSS = Twr.dO1_TSS;
dO2_TFA = Twr.dO2_TFA;
dO2_TSS = Twr.dO2_TSS;


nb = Bld.nb;
BldSec = Bld.BldSec;
O1_B1 = Bld.O1_B1;
O2_B1 = Bld.O2_B1;
O3_B1 = Bld.O3_B1;
W1_B1 = Bld.W1_B1;
W2_B1 = Bld.W2_B1;
W3_B1 = Bld.W3_B1;
s11_B1 = Bld.s11_B1;
s22_B1 = Bld.s22_B1;
s33_B1 = Bld.s33_B1;
s12_B1 = Bld.s12_B1;
s13_B1 = Bld.s13_B1;
s23_B1 = Bld.s23_B1;

Twist = Bld.Twist;
dO1_B1 = Bld.dO1_B1;
dO2_B1 = Bld.dO2_B1;
dO3_B1 = Bld.dO3_B1;
dW1_B1 = Bld.dW1_B1;
dW2_B1 = Bld.dW2_B1;
dW3_B1 = Bld.dW3_B1;

% Since the blades are similar the same properties can be copied
O1_B2 = O1_B1; O2_B2 = O2_B1; O3_B2 = O3_B1; W1_B2 = W1_B1; W2_B2 = W2_B1; W3_B2 = W3_B1;
s11_B2 = s11_B1; s22_B2 = s22_B1; s33_B2 = s33_B1; s12_B2 = s12_B1; s13_B2 = s13_B1; s23_B2 = s23_B1;

O1_B3 = O1_B1; O2_B3 = O2_B1; O3_B3 = O3_B1; W1_B3 = W1_B1; W2_B3 = W2_B1; W3_B3 = W3_B1;
s11_B3 = s11_B1; s22_B3 = s22_B1; s33_B3 = s33_B1; s12_B3 = s12_B1; s13_B3 = s13_B1; s23_B3 = s23_B1;

dO1_B2 = dO1_B1; dO2_B2 = dO2_B1; dO3_B2 = dO3_B1; dW1_B2 = dW1_B1; dW2_B2 = dW2_B1; dW3_B2 = dW3_B1;
dO1_B3 = dO1_B1; dO2_B3 = dO2_B1; dO3_B3 = dO3_B1; dW1_B3 = dW1_B1; dW2_B3 = dW2_B1; dW3_B3 = dW3_B1;


%% Generate results
for i = 1:length(t)
    BlPitch = Controls(i,2:4);
    %% Update state
    % Dofs
    Comq = zeros(DOFs.Avail*2,1);
    IndexActive = [DOFs.Active, DOFs.Active+DOFs.Avail];
    Comq(IndexActive) = q(i,:);
    q_Sg   = Comq(1);
    q_Sw   = Comq(2);
    q_Hv   = Comq(3);
    q_R    = Comq(4);
    q_P    = Comq(5);
    q_Y    = Comq(6);
    q_TFA1 = Comq(7);
    q_TSS1 = Comq(8);
    q_TFA2 = Comq(9);
    q_TSS2 = Comq(10);
    q_yaw  = Comq(11);
    q_GeAz = Comq(12);
    q_DrTr = Comq(13);
    q_B1F1 = Comq(14);
    q_B1E1 = Comq(15);
    q_B1F2 = Comq(16);
    q_B2F1 = Comq(17);
    q_B2E1 = Comq(18);
    q_B2F2 = Comq(19);
    q_B3F1 = Comq(20);
    q_B3E1 = Comq(21);
    q_B3F2 = Comq(22);

    qd_Sg   = Comq(23);
    qd_Sw   = Comq(24);
    qd_Hv   = Comq(25);
    qd_R    = Comq(26);
    qd_P    = Comq(27);
    qd_Y    = Comq(28);
    qd_TFA1 = Comq(29);
    qd_TSS1 = Comq(30);
    qd_TFA2 = Comq(31);
    qd_TSS2 = Comq(32);
    qd_yaw  = Comq(33);
    qd_GeAz = Comq(34);
    qd_DrTr = Comq(35);
    qd_B1F1 = Comq(36);
    qd_B1E1 = Comq(37);
    qd_B1F2 = Comq(38);
    qd_B2F1 = Comq(39);
    qd_B2E1 = Comq(40);
    qd_B2F2 = Comq(41);
    qd_B3F1 = Comq(42);
    qd_B3E1 = Comq(43);
    qd_B3F2 = Comq(44);
    % end of Update state

    %% Coordinate Systems

    Z = eye(3);                                                                 % Inertial/Earth coordinate system

    A = Transform1(FASTTransMat(q_R, q_Y, -q_P),Z);                                             % Tower coordinate system

    theta_FA = -dO1_TFA*q_TFA1 - dO2_TFA*q_TFA2;
    theta_SS =  dO1_TSS*q_TSS1 + dO2_TSS*q_TSS2;

    B = Transform1(FASTTransMat(theta_SS, 0, theta_FA),A);                                      % Tower top/ Base plate coordinate system

    D = Transform1([cos(q_yaw), 0, -sin(q_yaw); 0 1 0; sin(q_yaw) 0 cos(q_yaw)],B);         % Nacelle coordinate system
    ShftSkew = 0;

    C = Transform1([cos(ShftSkew)*cos(ShftTilt), sin(ShftTilt), -sin(ShftSkew)*cos(ShftTilt);                            % Shaft coordinate system
        -cos(ShftSkew)*sin(ShftTilt), cos(ShftTilt),  sin(ShftSkew)*sin(ShftTilt);
         sin(ShftSkew),                     0,         cos(ShftSkew)],D);

    E = Transform1([1, 0, 0; 0, cos(q_DrTr+q_GeAz), sin(q_DrTr+q_GeAz);                   % Azimuth coordinate system
          0, -sin(q_DrTr+q_GeAz), cos(q_DrTr+q_GeAz)],C);

    G1 = Transform1([1, 0, 0; 0, cos(-pi/2), sin(-pi/2);                                    % Coordinate system for blade 1
          0, -sin(-pi/2), cos(-pi/2)],E);

    G2 = Transform1([1, 0, 0; 0, cos(-pi/2 + 2*pi/3), sin(-pi/2 + 2*pi/3);                  % Coordinate system for blade 1
          0, -sin(-pi/2 + 2*pi/3), cos(-pi/2 + 2*pi/3)],E);

    G3 = Transform1([1, 0, 0; 0, cos(-pi/2 + 4*pi/3), sin(-pi/2 + 4*pi/3);                  % Coordinate system for blade 1
          0, -sin(-pi/2 + 4*pi/3), cos(-pi/2 + 4*pi/3)],E);

    I1 = Transform1([cos(PreCone(1)), 0, -sin(PreCone(1));
                0,         1,       0;                                          % Coned coordinate system for blade 1
          sin(PreCone(1)), 0, cos(PreCone(1))],G1);

    I2 = Transform1([cos(PreCone(2)), 0, -sin(PreCone(2));
                0,         1,       0;                                          % Coned coordinate system for blade 1
          sin(PreCone(2)), 0, cos(PreCone(2))],G2);

    I3 = Transform1([cos(PreCone(3)), 0, -sin(PreCone(3));
                0,         1,       0;                                          % Coned coordinate system for blade 1
          sin(PreCone(3)), 0, cos(PreCone(3))],G3);

    J1 = Transform1([cos(BlPitch(1)), -sin(BlPitch(1)), 0;                                 % Pitched coordinate system for blade 1
          sin(BlPitch(1)),  cos(BlPitch(1)), 0;
                0,                 0,        1],I1);

    J2 = Transform1([cos(BlPitch(2)), -sin(BlPitch(2)), 0;                                 % Pitched coordinate system for blade 1
          sin(BlPitch(2)),  cos(BlPitch(2)), 0;
                0,                 0,        1],I2);

    J3 = Transform1([cos(BlPitch(3)), -sin(BlPitch(3)), 0;                                 % Pitched coordinate system for blade 1
          sin(BlPitch(3)),  cos(BlPitch(3)), 0;
                0,                 0,        1],I3);

    % Twist = reshape(Twist,1,1,nb_Load);
    % 
    % LB1 = Transform2([cos(Twist), -sin(Twist), zeros(1,1,nb_Load);
    %                    sin(Twist),  cos(Twist), zeros(1,1,nb_Load);
    %                zeros(1,1,nb_Load), zeros(1,1,nb_Load), ones(1,1,nb_Load)],J1);
    %            
    % LB2 = Transform2([cos(Twist), -sin(Twist), zeros(1,1,nb_Load);
    %                    sin(Twist),  cos(Twist), zeros(1,1,nb_Load);
    %                zeros(1,1,nb_Load), zeros(1,1,nb_Load), ones(1,1,nb_Load)],J2);
    %            
    % LB3 = Transform2([cos(Twist), -sin(Twist), zeros(1,1,nb_Load);
    %                    sin(Twist),  cos(Twist), zeros(1,1,nb_Load);
    %                zeros(1,1,nb_Load), zeros(1,1,nb_Load), ones(1,1,nb_Load)],J3);
    %    
    % % Blade Element Fixed Coordinate System
    % thetaB1_IP  = -(dW1_B1*q_B1F1 + dW2_B1*q_B1F2 + dW3_B1*q_B1E1);
    % thetaB1_Oop =  (dO1_B1*q_B1F1 + dO2_B1*q_B1F2 + dO3_B1*q_B1E1);
    % 
    % thetaB1_x = cos(Twist).*thetaB1_IP - sin(Twist).*thetaB1_Oop;
    % thetaB1_y = sin(Twist).*thetaB1_IP + cos(Twist).*thetaB1_Oop;
    % 
    % thetaB2_IP  = -(dW1_B2*q_B2F1 + dW2_B2*q_B2F2 + dW3_B2*q_B2E1);
    % thetaB2_Oop =  (dO1_B2*q_B2F1 + dO2_B2*q_B2F2 + dO3_B2*q_B2E1);
    % 
    % thetaB2_x = cos(Twist).*thetaB2_IP - sin(Twist).*thetaB2_Oop;
    % thetaB2_y = sin(Twist).*thetaB2_IP + cos(Twist).*thetaB2_Oop;
    % 
    % thetaB3_IP  = -(dW1_B3*q_B3F1 + dW2_B3*q_B3F2 + dW3_B3*q_B3E1);
    % thetaB3_Oop =  (dO1_B3*q_B3F1 + dO2_B3*q_B3F2 + dO3_B3*q_B3E1);
    % 
    % thetaB3_x = cos(Twist).*thetaB3_IP - sin(Twist).*thetaB3_Oop;
    % thetaB3_y = sin(Twist).*thetaB3_IP + cos(Twist).*thetaB3_Oop;
    % 
    % NB1 = Transform3(FASTTransMat(thetaB1_x, thetaB1_y, zeros(1,1,nb_Load)),LB1);            % TransMat must return a 3D matrix for vectors thetaB1_x and thetaB1_y
    % NB2 = Transform3(FASTTransMat(thetaB2_x, thetaB2_y, zeros(1,1,nb_Load)),LB2);
    % NB3 = Transform3(FASTTransMat(thetaB3_x, thetaB3_y, zeros(1,1,nb_Load)),LB3);
    % 
    % % Blade Element Fixed Coordinate System for Returning Aerodynamic Loads
    % 
    % MB1 = Transform3([cos(Twist + BlPitch(1)),  sin(Twist + BlPitch(1)), zeros(1,1,nb_Load) ;
    %                   -sin(Twist + BlPitch(1)),  cos(Twist + BlPitch(1)), zeros(1,1,nb_Load);
    %                    zeros(1,1,nb_Load),        zeros(1,1,nb_Load),     ones(1,1,nb_Load)],NB1);
    %        
    % MB2 = Transform3([cos(Twist + BlPitch(2)),  sin(Twist + BlPitch(2)), zeros(1,1,nb_Load) ;
    %                   -sin(Twist + BlPitch(2)),  cos(Twist + BlPitch(2)), zeros(1,1,nb_Load);
    %                    zeros(1,1,nb_Load),        zeros(1,1,nb_Load),     ones(1,1,nb_Load)],NB2);
    %        
    % MB3 = Transform3([cos(Twist + BlPitch(3)),  sin(Twist + BlPitch(3)), zeros(1,1,nb_Load) ;
    %                   -sin(Twist + BlPitch(3)),  cos(Twist + BlPitch(3)), zeros(1,1,nb_Load);
    %                    zeros(1,1,nb_Load),        zeros(1,1,nb_Load),     ones(1,1,nb_Load)],NB3);

    %% Position Vectors

    rZ  = q_Sg*Z(1,:) + q_Hv*Z(2,:) - q_Sw*Z(3,:);                               % platform reference
    rZY = (PtfmRefzt-PtfmCMzt)*A(2,:);                                               % platform reference to platform C.M

    rZT_1 = O1_TFA*q_TFA1 + O2_TFA*q_TFA2;                                     % O1_TF, O2_TF, s11_TFA, s12_TFA etc are vectors
    rZT_2 = reshape(TwrSec,1,1,nt) + PtfmRefzt + TowerBsHt - 0.5*(s11_TFA*q_TFA1^2 + s22_TFA ...        
               *q_TFA2^2 +2*s12_TFA*q_TFA1*q_TFA2 + s11_TSS*q_TSS1^2 + s22_TSS*q_TSS2^2 ...
                + 2*s12_TSS*q_TSS1*q_TSS2);        
    rZT_3 = O1_TSS*q_TSS1 + O2_TSS*q_TSS2;
    rZT = rZT_1.*A(1,:) + rZT_2.*A(2,:) + rZT_3.*A(3,:);                         % rZT is a 3D matrix

    rZO_1 = q_TFA1 + q_TFA2;                                                     % S11_TFA, S22_TF etc are scalars
    rZO_2 = PtfmRefzt + TwrHt - 0.5*(S11_TFA*q_TFA1^2 + S22_TFA ...
               *q_TFA2^2 +2*S12_TFA*q_TFA1*q_TFA2 + S11_TSS*q_TSS1^2 + S22_TSS*q_TSS2^2 ...
                + 2*S12_TSS*q_TSS1*q_TSS2);
    rZO_3 = q_TSS1 + q_TSS2;
    rZO   = rZO_1*A(1,:) + rZO_2*A(2,:) + rZO_3*A(3,:); 

%     rOU = NacCMxn*D(1,:) + NacCMzn*D(2,:) - NacCMyn*D(3,:);
    % rOIMU = NcIMUxn*D(1,:) + NcIMUzn*D(2,:) - NcIMUyn*D(3,:);
    Yaw2Shft =0;
    rOQ = OverHang*C(1,:) + Twr2Shft*D(2,:) - Yaw2Shft*D(3,:);
%     rQC = HubCM*C(1,:);

    rQS1_1 = O1_B1*q_B1F1 + O2_B1*q_B1F2 + O3_B1*q_B1E1;
    rQS1_3 = reshape(BldSec,1,1,nb) + HubRad - 0.5*(s11_B1*q_B1F1^2 + s22_B1*q_B1F2^2 + s33_B1*q_B1E1^2 ...
           + 2*s12_B1*q_B1F1*q_B1F2 + 2*s23_B1*q_B1F2*q_B1E1 + 2*s13_B1*q_B1F1*q_B1E1);
    rQS1_2 = W1_B1*q_B1F1 + W2_B1*q_B1F2 + W3_B1*q_B1E1;

    rQS1   = rQS1_1.*J1(1,:) + rQS1_2.*J1(2,:) + rQS1_3.*J1(3,:);

    rQS2_1 = O1_B2*q_B2F1 + O2_B2*q_B2F2 + O3_B2*q_B2E1;
    rQS2_3 = reshape(BldSec,1,1,nb) + HubRad - 0.5*(s11_B2*q_B2F1^2 + s22_B2*q_B2F2^2 + s33_B2*q_B2E1^2 ...
           + 2*s12_B2*q_B2F1*q_B2F2 + 2*s23_B2*q_B2F2*q_B2E1 + 2*s13_B2*q_B2F1*q_B2E1);
    rQS2_2 = W1_B2*q_B2F1 + W2_B2*q_B2F2 + W3_B2*q_B2E1;

    rQS2 = rQS2_1.*J2(1,:) + rQS2_2.*J2(2,:) + rQS2_3.*J2(3,:);

    rQS3_1 = O1_B3*q_B3F1 + O2_B3*q_B3F2 + O3_B3*q_B3E1;
    rQS3_3 = reshape(BldSec,1,1,nb) + HubRad - 0.5*(s11_B3*q_B3F1^2 + s22_B3*q_B3F2^2 + s33_B3*q_B3E1^2 ...
           + 2*s12_B3*q_B3F1*q_B3F2 + 2*s23_B3*q_B3F2*q_B3E1 + 2*s13_B3*q_B3F1*q_B3E1);
    rQS3_2 = W1_B3*q_B3F1 + W2_B3*q_B3F2 + W3_B3*q_B3E1;

    rQS3 = rQS3_1.*J3(1,:) + rQS3_2.*J3(2,:) + rQS3_3.*J3(3,:);

    %% New vectors
    rZZT = rZ + rZT(1,:,5);
    rZS1 = rZ + rZO + rOQ + rQS1(1,:,end);
    
    vectordist = rZS1 - rZZT;
    
    %% Outputs
    OopDefl1(i) = dot(rQS1(1,:,end) - TipRad*J1(3,:),I1(1,:));
    IPDefl1(i)  = dot(rQS1(1,:,end) - TipRad*J1(3,:),I1(2,:));
    TipDzb1(i)  = dot(rQS1(1,:,end) - TipRad*J1(3,:),I1(3,:));
    
    TTDspFA(i)  = dot(rZO - (TwrHt+PtfmRefzt)*A(2,:),A(1,:));
    TTDspSS(i)  = -dot(rZO - (TwrHt+PtfmRefzt)*A(2,:),A(3,:));
    
    LSSTipV(i)  = 30/pi*(qd_DrTr+qd_GeAz);
    NacYaw(i)   = 180/pi*q_yaw;
    
    % Tower Blade 1 Clearance
    rQB1FlexLen = rQS1(1,:,end);
    rOS1BlFlexLen = rOQ + rQB1FlexLen;
    if dot(rOS1BlFlexLen, D(2,:)) > 0
        TwrClrnB1(i) = sqrt( dot(rOS1BlFlexLen, D(1,:))^2 + dot(rOS1BlFlexLen, D(2,:))^2 + dot(rOS1BlFlexLen, D(3,:))^2 );
    else
        TwrClrnB1(i) = sqrt( dot(rOS1BlFlexLen, D(1,:))^2 + dot(rOS1BlFlexLen, D(3,:))^2 );
    end
    
    % Tower Blade 2 Clearance
    rQB2FlexLen = rQS2(1,:,end);
    rOS2BlFlexLen = rOQ + rQB2FlexLen;   
    if dot(rOS2BlFlexLen, D(2,:)) > 0
        TwrClrnB2(i) = sqrt( dot(rOS2BlFlexLen, D(1,:))^2 + dot(rOS2BlFlexLen, D(2,:))^2 + dot(rOS2BlFlexLen, D(3,:))^2 );
    else
        TwrClrnB2(i) = sqrt( dot(rOS2BlFlexLen, D(1,:))^2 + dot(rOS2BlFlexLen, D(3,:))^2 );
    end
    
    % Tower Blade 3 Clearance
    rQB3FlexLen = rQS3(1,:,end);  
    rOS3BlFlexLen = rOQ + rQB3FlexLen;    
    if dot(rOS3BlFlexLen, D(2,:)) > 0
        TwrClrnB3(i) = sqrt( dot(rOS3BlFlexLen, D(1,:))^2 + dot(rOS3BlFlexLen, D(2,:))^2 + dot(rOS3BlFlexLen, D(3,:))^2 );
    else
        TwrClrnB3(i) = sqrt( dot(rOS3BlFlexLen, D(1,:))^2 + dot(rOS3BlFlexLen, D(3,:))^2 );
    end
    
    Surge(i) = q_Sg;
    Sway(i)  = q_Sw;
    Heave(i) = q_Hv;
    Roll(i)  = q_R;
    Pitch(i) = q_P;
    Yaw(i)   = q_Y;
    
    BldTwrHrClr(i) = vectordist(1);
end
end

