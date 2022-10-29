function [Z, A, B, C, E, J1, J2, J3, MB1, MB2, MB3] = Coordinate_systems(q_Nom, BlPitch, ElastoDyn, Twr, Bld)
PreCone = ElastoDyn.PreCone;
ShftTilt = ElastoDyn.ShftTilt;
dO1_B1 = Bld.dO1_B1;
dO2_B1 = Bld.dO2_B1;
dO3_B1 = Bld.dO3_B1;
dW1_B1 = Bld.dW1_B1;
dW2_B1 = Bld.dW2_B1;
dW3_B1 = Bld.dW3_B1;
Twist = Bld.Twist;
dO1_TFA = Twr.dO1_TFA;
dO1_TSS = Twr.dO1_TSS;
dO2_TFA = Twr.dO2_TFA;
dO2_TSS = Twr.dO2_TSS;

dO1_B2 = dO1_B1; dO2_B2 = dO2_B1; dO3_B2 = dO3_B1; dW1_B2 = dW1_B1; dW2_B2 = dW2_B1; dW3_B2 = dW3_B1;
dO1_B3 = dO1_B1; dO2_B3 = dO2_B1; dO3_B3 = dO3_B1; dW1_B3 = dW1_B1; dW2_B3 = dW2_B1; dW3_B3 = dW3_B1;

%% Dofs
% q_Sg   = q_Nom(1);
% q_Sw   = q_Nom(2);
% q_Hv   = q_Nom(3);
q_R    = q_Nom(4);
q_P    = q_Nom(5);
q_Y    = q_Nom(6);
q_TFA1 = q_Nom(7);
q_TSS1 = q_Nom(8);
q_TFA2 = q_Nom(9);
q_TSS2 = q_Nom(10);
q_yaw  = q_Nom(11);
q_GeAz = q_Nom(12);
q_DrTr = q_Nom(13);
q_B1F1 = q_Nom(14);
q_B1E1 = q_Nom(15);
q_B1F2 = q_Nom(16);
q_B2F1 = q_Nom(17);
q_B2E1 = q_Nom(18);
q_B2F2 = q_Nom(19);
q_B3F1 = q_Nom(20);
q_B3E1 = q_Nom(21);
q_B3F2 = q_Nom(22);
 
% qd_Sg   = q_Nom(23);
% qd_Sw   = q_Nom(24);
% qd_Hv   = q_Nom(25);
% qd_R    = q_Nom(26);
% qd_P    = q_Nom(27);
% qd_Y    = q_Nom(28);
% qd_TFA1 = q_Nom(29);
% qd_TSS1 = q_Nom(30);
% qd_TFA2 = q_Nom(31);
% qd_TSS2 = q_Nom(32);
% qd_yaw  = q_Nom(33);
% qd_GeAz = q_Nom(34);
% qd_DrTr = q_Nom(35);
% qd_B1F1 = q_Nom(36);
% qd_B1E1 = q_Nom(37);
% qd_B1F2 = q_Nom(38);
% qd_B2F1 = q_Nom(39);
% qd_B2E1 = q_Nom(40);
% qd_B2F2 = q_Nom(41);
% qd_B3F1 = q_Nom(42);
% qd_B3E1 = q_Nom(43);
% qd_B3F2 = q_Nom(44);
% end of Update state

Z = eye(3);                                                                 % Inertial/Earth coordinate system

A = Transform1(FASTTransMat(q_R, q_Y, -q_P), Z);                                             % Tower coordinate system

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
nb = Bld.nb;
Twist = reshape(Twist,1,1,nb);

LB1 = Transform2([cos(Twist), -sin(Twist), zeros(1,1,nb);
                   sin(Twist),  cos(Twist), zeros(1,1,nb);
               zeros(1,1,nb), zeros(1,1,nb), ones(1,1,nb)],J1);
           
LB2 = Transform2([cos(Twist), -sin(Twist), zeros(1,1,nb);
                   sin(Twist),  cos(Twist), zeros(1,1,nb);
               zeros(1,1,nb), zeros(1,1,nb), ones(1,1,nb)],J2);
           
LB3 = Transform2([cos(Twist), -sin(Twist), zeros(1,1,nb);
                   sin(Twist),  cos(Twist), zeros(1,1,nb);
               zeros(1,1,nb), zeros(1,1,nb), ones(1,1,nb)],J3);
   
% Blade Element Fixed Coordinate System

thetaB1_IP  = -(dW1_B1*q_B1F1 + dW2_B1*q_B1F2 + dW3_B1*q_B1E1);
thetaB1_Oop =  (dO1_B1*q_B1F1 + dO2_B1*q_B1F2 + dO3_B1*q_B1E1);

thetaB1_x = cos(Twist).*thetaB1_IP - sin(Twist).*thetaB1_Oop;
thetaB1_y = sin(Twist).*thetaB1_IP + cos(Twist).*thetaB1_Oop;

thetaB2_IP  = -(dW1_B2*q_B2F1 + dW2_B2*q_B2F2 + dW3_B2*q_B2E1);
thetaB2_Oop =  (dO1_B2*q_B2F1 + dO2_B2*q_B2F2 + dO3_B2*q_B2E1);

thetaB2_x = cos(Twist).*thetaB2_IP - sin(Twist).*thetaB2_Oop;
thetaB2_y = sin(Twist).*thetaB2_IP + cos(Twist).*thetaB2_Oop;

thetaB3_IP  = -(dW1_B3*q_B3F1 + dW2_B3*q_B3F2 + dW3_B3*q_B3E1);
thetaB3_Oop =  (dO1_B3*q_B3F1 + dO2_B3*q_B3F2 + dO3_B3*q_B3E1);

thetaB3_x = cos(Twist).*thetaB3_IP - sin(Twist).*thetaB3_Oop;
thetaB3_y = sin(Twist).*thetaB3_IP + cos(Twist).*thetaB3_Oop;

NB1 = Transform3(FASTTransMat(thetaB1_x, thetaB1_y, zeros(1,1,nb)),LB1);    % FASTTransMat must return a 3D matrix for vectors thetaB1_x and thetaB1_y
NB2 = Transform3(FASTTransMat(thetaB2_x, thetaB2_y, zeros(1,1,nb)),LB2);
NB3 = Transform3(FASTTransMat(thetaB3_x, thetaB3_y, zeros(1,1,nb)),LB3);

% Blade Element Fixed Coordinate System for Returning Aerodynamic Loads

MB1 = Transform3([cos(Twist + BlPitch(1)),  sin(Twist + BlPitch(1)), zeros(1,1,nb) ;
                  -sin(Twist + BlPitch(1)),  cos(Twist + BlPitch(1)), zeros(1,1,nb);
                   zeros(1,1,nb),        zeros(1,1,nb),     ones(1,1,nb)],NB1);
       
MB2 = Transform3([cos(Twist + BlPitch(2)),  sin(Twist + BlPitch(2)), zeros(1,1,nb) ;
                  -sin(Twist + BlPitch(2)),  cos(Twist + BlPitch(2)), zeros(1,1,nb);
                   zeros(1,1,nb),        zeros(1,1,nb),     ones(1,1,nb)],NB2);
       
MB3 = Transform3([cos(Twist + BlPitch(3)),  sin(Twist + BlPitch(3)), zeros(1,1,nb) ;
                  -sin(Twist + BlPitch(3)),  cos(Twist + BlPitch(3)), zeros(1,1,nb);
                   zeros(1,1,nb),        zeros(1,1,nb),     ones(1,1,nb)],NB3);
               
end

