function [IM_nom, f_nom, Controls] = NominalSystemMatrix(q_Nom, Controls, ElastoDyn, Airfoils, Twr, Bld, Platform, WindNom, mooring_load, f_Morison)
%#Codegen
g = 9.80665;
% Dofs
q_Sg   = q_Nom(1);
q_Sw   = q_Nom(2);
q_Hv   = q_Nom(3);
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
 
qd_Sg   = q_Nom(23);
qd_Sw   = q_Nom(24);
qd_Hv   = q_Nom(25);
qd_R    = q_Nom(26);
qd_P    = q_Nom(27);
qd_Y    = q_Nom(28);
qd_TFA1 = q_Nom(29);
qd_TSS1 = q_Nom(30);
qd_TFA2 = q_Nom(31);
qd_TSS2 = q_Nom(32);
qd_yaw  = q_Nom(33);
qd_GeAz = q_Nom(34);
qd_DrTr = q_Nom(35);
qd_B1F1 = q_Nom(36);
qd_B1E1 = q_Nom(37);
qd_B1F2 = q_Nom(38);
qd_B2F1 = q_Nom(39);
qd_B2E1 = q_Nom(40);
qd_B2F2 = q_Nom(41);
qd_B3F1 = q_Nom(42);
qd_B3E1 = q_Nom(43);
qd_B3F2 = q_Nom(44);
% end of Update state

%% Read Data from Structure to Doubles
TipRad = ElastoDyn.TipRad;
HubRad = ElastoDyn.HubRad;
HubCM = ElastoDyn.HubCM;
OverHang = ElastoDyn.OverHang;
NacCMxn = ElastoDyn.NacCMxn;
NacCMyn = ElastoDyn.NacCMyn;
NacCMzn = ElastoDyn.NacCMzn;
Twr2Shft = ElastoDyn.Twr2Shft;
TwrHt = ElastoDyn.TwrHt;
TowerBsHt = ElastoDyn.TowerBsHt;  
PtfmCMzt = ElastoDyn.PtfmCMzt;   
PtfmRefzt = ElastoDyn.PtfmRefzt;   

mH = ElastoDyn.HubMass;
HubIner = ElastoDyn.HubIner;
GenIner = ElastoDyn.GenIner;
mN = ElastoDyn.NacMass;
NacYIner = ElastoDyn.NacYIner;
YawBrMass = ElastoDyn.YawBrMass;
mX = ElastoDyn.PtfmMass;
PtfmRIner = ElastoDyn.PtfmRIner;
PtfmPIner = ElastoDyn.PtfmPIner;
PtfmYIner = ElastoDyn.PtfmYIner;

GenDir  = 1;
GBRatio = ElastoDyn.GBRatio;
%GBoxEff = ElastoDyn.GBoxEff;
DTTorSpr = ElastoDyn.DTTorSpr;
DTTorDmp = ElastoDyn.DTTorDmp;

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
mT = Twr.mT;
k11_TFA = Twr.k11_TFA;
k12_TFA = Twr.k12_TFA;
k11_TSS = Twr.k11_TSS;
k12_TSS = Twr.k12_TSS;
k21_TFA = Twr.k21_TFA;
k22_TFA = Twr.k22_TFA;
k21_TSS = Twr.k21_TSS;
k22_TSS = Twr.k22_TSS;
f1_TFA = Twr.f1_TFA;
f2_TFA = Twr.f2_TFA;
f1_TSS = Twr.f1_TSS;
f2_TSS = Twr.f2_TSS;
zeta1_TFA = Twr.zeta1_TFA;
zeta2_TFA = Twr.zeta2_TFA;
zeta1_TSS = Twr.zeta1_TSS;
zeta2_TSS = Twr.zeta2_TSS;

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
mB1 = Bld.mB1;
k11_B1F = Bld.k11_B1F;
k12_B1F = Bld.k12_B1F;
k21_B1F = k12_B1F;
k22_B1F = Bld.k22_B1F;
k11_B1E = Bld.k11_B1E;
f1_B1F = Bld.f1_B1F;
f2_B1F = Bld.f2_B1F;
f1_B1E = Bld.f1_B1E;
zeta1_B1F = Bld.zeta1_B1F;
zeta2_B1F = Bld.zeta2_B1F;
zeta1_B1E = Bld.zeta1_B1E;

YawNeut = 0;
YawSpr  = 9.02832E+09;
YawDamp = 1.916E+07;

% BlSpn = Load.BlSpn;
% nb = Load.nb;
Cord = Bld.Cord;
FoilNo = Bld.FoilNo;
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

k11_B2F = k11_B1F; k12_B2F = k12_B1F; k22_B2F = k22_B1F; k11_B2E = k11_B1E; k21_B2F = k12_B2F; f1_B2F = f1_B1F; f2_B2F = f2_B1F; f1_B2E = f1_B1E;
k11_B3F = k11_B1F; k12_B3F = k12_B1F; k22_B3F = k22_B1F; k11_B3E = k11_B1E; k21_B3F = k12_B3F; f1_B3F = f1_B1F; f2_B3F = f2_B1F; f1_B3E = f1_B1E;

zeta1_B2F = zeta1_B1F; zeta2_B2F = zeta2_B1F; zeta1_B2E = zeta1_B1E; zeta1_B3F = zeta1_B1F; zeta2_B3F = zeta2_B1F; zeta1_B3E = zeta1_B1E;

mB2 = mB1; mB3 = mB1;

% TMDzn    = TMD.TMDzn;

% end of reading

%% Coordinate Systems
BlPitch = Controls(2:4);
GenTrq  = Controls(1);
[Z, A, D, C, E, J1, J2, J3, ~, ~, ~, MB1, MB2, MB3] = Coordinate_systems(q_Nom, BlPitch, ElastoDyn, Twr, Bld);

%% Position vectors
% Position Vectors

% rZ  = q_Sg*Z(1,:) + q_Hv*Z(2,:) - q_Sw*Z(3,:);                               % platform reference
rZY = (PtfmCMzt - PtfmRefzt)*A(2,:);                                               % platform reference to platform C.M

rZT_1 = O1_TFA*q_TFA1 + O2_TFA*q_TFA2;                                     % O1_TF, O2_TF, s11_TFA, s12_TFA etc are vectors
rZT_2 = reshape(TwrSec,1,1,nt) + PtfmRefzt + TowerBsHt - 0.5*(s11_TFA*q_TFA1^2 + s22_TFA ...        
           *q_TFA2^2 +2*s12_TFA*q_TFA1*q_TFA2 + s11_TSS*q_TSS1^2 + s22_TSS*q_TSS2^2 ...
            + 2*s12_TSS*q_TSS1*q_TSS2);        
rZT_3 = O1_TSS*q_TSS1 + O2_TSS*q_TSS2;
rZT = coprod(rZT_1,A(1,:)) + coprod(rZT_2,A(2,:)) + coprod(rZT_3,A(3,:));                         % rZT is a 3D matrix
rZT0 = rZT(:,:,1);

rZO_1 = q_TFA1 + q_TFA2;                                                     % S11_TFA, S22_TF etc are scalars
rZO_2 = PtfmRefzt + TwrHt - 0.5*(S11_TFA*q_TFA1^2 + S22_TFA ...
           *q_TFA2^2 +2*S12_TFA*q_TFA1*q_TFA2 + S11_TSS*q_TSS1^2 + S22_TSS*q_TSS2^2 ...
            + 2*S12_TSS*q_TSS1*q_TSS2);
rZO_3 = q_TSS1 + q_TSS2;
rZO   = rZO_1*A(1,:) + rZO_2*A(2,:) + rZO_3*A(3,:); 

rOU = NacCMxn*D(1,:) + NacCMzn*D(2,:) - NacCMyn*D(3,:);

% rOIMU = NcIMUxn*D(1,:) + NcIMUzn*D(2,:) - NcIMUyn*D(3,:);
Yaw2Shft =0;
rOQ = OverHang*C(1,:) + Twr2Shft*D(2,:) - Yaw2Shft*D(3,:);

rQC = HubCM*E(1,:);

rQS1_1 = O1_B1*q_B1F1 + O2_B1*q_B1F2 + O3_B1*q_B1E1;
rQS1_3 = reshape(BldSec,1,1,nb) + HubRad - 0.5*(s11_B1*q_B1F1^2 + s22_B1*q_B1F2^2 + s33_B1*q_B1E1^2 ...
       + 2*s12_B1*q_B1F1*q_B1F2 + 2*s23_B1*q_B1F2*q_B1E1 + 2*s13_B1*q_B1F1*q_B1E1);
rQS1_2 = W1_B1*q_B1F1 + W2_B1*q_B1F2 + W3_B1*q_B1E1;

rQS1 = coprod(rQS1_1,J1(1,:)) + coprod(rQS1_2,J1(2,:)) + coprod(rQS1_3,J1(3,:));
rQS10 = rQS1(:,:,1);

rQS2_1 = O1_B2*q_B2F1 + O2_B2*q_B2F2 + O3_B2*q_B2E1;
rQS2_3 = reshape(BldSec,1,1,nb) + HubRad - 0.5*(s11_B2*q_B2F1^2 + s22_B2*q_B2F2^2 + s33_B2*q_B2E1^2 ...
       + 2*s12_B2*q_B2F1*q_B2F2 + 2*s23_B2*q_B2F2*q_B2E1 + 2*s13_B2*q_B2F1*q_B2E1);
rQS2_2 = W1_B2*q_B2F1 + W2_B2*q_B2F2 + W3_B2*q_B2E1;
rQS2 = coprod(rQS2_1,J2(1,:)) + coprod(rQS2_2,J2(2,:)) + coprod(rQS2_3,J2(3,:));
rQS20 = rQS2(:,:,1);

rQS3_1 = O1_B3*q_B3F1 + O2_B3*q_B3F2 + O3_B3*q_B3E1;
rQS3_3 = reshape(BldSec,1,1,nb) + HubRad - 0.5*(s11_B3*q_B3F1^2 + s22_B3*q_B3F2^2 + s33_B3*q_B3E1^2 ...
       + 2*s12_B3*q_B3F1*q_B3F2 + 2*s23_B3*q_B3F2*q_B3E1 + 2*s13_B3*q_B3F1*q_B3E1);
rQS3_2 = W1_B3*q_B3F1 + W2_B3*q_B3F2 + W3_B3*q_B3E1;

rQS3 = coprod(rQS3_1,J3(1,:)) + coprod(rQS3_2,J3(2,:)) + coprod(rQS3_3,J3(3,:));
rQS30 = rQS3(:,:,1);
% end of Position Vectors

%% Angular velocities

EwX = qd_R*Z(1,:) + qd_Y*Z(2,:) - qd_P*Z(3,:);
EwB = EwX + (dO1_TSS*qd_TSS1 + dO2_TSS*qd_TSS2)*A(1,:) - (dO1_TFA*qd_TFA1 + dO2_TFA*qd_TFA2)*A(3,:);
EwN = EwB + qd_yaw*D(2,:);
EwL = EwN + qd_DrTr*C(1,:) + qd_GeAz*C(1,:);
EwG = EwN + GenDir*GBRatio*qd_GeAz*C(1,:);
% end of Angular velocities

%% Velocity vectors

EvZ   = qd_Sg*Z(1,:) + qd_Hv*Z(2,:) - qd_Sw*Z(3,:);

XvT_1 = O1_TFA*qd_TFA1 + O2_TFA*qd_TFA2;
XvT_2 = -(s11_TFA*q_TFA1*qd_TFA1 + s22_TFA*q_TFA2*qd_TFA2 + s12_TFA*(qd_TFA1*q_TFA2+q_TFA1*qd_TFA2) ...
        + s11_TSS*q_TSS1*qd_TSS1 + s22_TSS*q_TSS2*qd_TSS2 + s12_TSS*(qd_TSS1*q_TSS2+q_TSS1*qd_TSS2));
XvT_3 = O1_TSS*qd_TSS1 + O2_TSS*qd_TSS2;
XvT = coprod(XvT_1,A(1,:)) + coprod(XvT_2,A(2,:)) + coprod(XvT_3,A(3,:));

XvO_1 = qd_TFA1 + qd_TFA2;
XvO_2 = -(S11_TFA*q_TFA1*qd_TFA1 + S22_TFA*q_TFA2*qd_TFA2 + S12_TFA*(qd_TFA1*q_TFA2+q_TFA1*qd_TFA2) ...
        + S11_TSS*q_TSS1*qd_TSS1 + S22_TSS*q_TSS2*qd_TSS2 + S12_TSS*(qd_TSS1*q_TSS2+q_TSS1*qd_TSS2));
XvO_3 = qd_TSS1 + qd_TSS2;
XvO = XvO_1*A(1,:) + XvO_2*A(2,:) + XvO_3*A(3,:);

EvO    = XvO + EvZ + cross(EwX,rZO);
EvQ    = EvO + cross(EwN,rOQ);
HvS1_1 = O1_B1*qd_B1F1 + O2_B1*qd_B1F2 + O3_B1*qd_B1E1;
HvS1_3 = -(s11_B1*q_B1F1*qd_B1F1 + s22_B1*q_B1F2*qd_B1F2 + s33_B1*q_B1E1*qd_B1E1 ...
        + s12_B1*(qd_B1F1*q_B1F2 + q_B1F1*qd_B1F2) + s23_B1*(qd_B1F2*q_B1E1 + q_B1F2*qd_B1E1) ...
        + s13_B1*(qd_B1F1*q_B1E1 + q_B1F1*qd_B1E1));
HvS1_2 = W1_B1*qd_B1F1 + W2_B1*qd_B1F2 + W3_B1*qd_B1E1;
HvS1 = coprod(HvS1_1,J1(1,:)) + coprod(HvS1_2,J1(2,:)) + coprod(HvS1_3,J1(3,:));

EvS1   =  repmat(EvQ,1,1,nb) + HvS1 + cross(repmat(EwL,1,1,nb),rQS1,2); 

HvS2_1 = O1_B2*qd_B2F1 + O2_B2*qd_B2F2 + O3_B2*qd_B2E1;
HvS2_3 = -(s11_B2*q_B2F1*qd_B2F1 + s22_B2*q_B2F2*qd_B2F2 + s33_B2*q_B2E1*qd_B2E1 ...
        + s12_B2*(qd_B2F1*q_B2F2 + q_B2F1*qd_B2F2) + s23_B2*(qd_B2F2*q_B2E1 + q_B2F2*qd_B2E1) ...
        + s13_B2*(qd_B2F1*q_B2E1 + q_B2F1*qd_B2E1));
HvS2_2 = W1_B2*qd_B2F1 + W2_B2*qd_B2F2 + W3_B2*qd_B2E1;
HvS2 = coprod(HvS2_1,J2(1,:)) + coprod(HvS2_2,J2(2,:)) + coprod(HvS2_3,J2(3,:));

EvS2   =  repmat(EvQ,1,1,nb) + HvS2 + cross(repmat(EwL,1,1,nb),rQS2,2);

HvS3_1 = O1_B3*qd_B3F1 + O2_B3*qd_B3F2 + O3_B3*qd_B3E1;
HvS3_3 = -(s11_B3*q_B3F1*qd_B3F1 + s22_B3*q_B3F2*qd_B3F2 + s33_B3*q_B3E1*qd_B3E1 ...
        + s12_B3*(qd_B3F1*q_B3F2 + q_B3F1*qd_B3F2) + s23_B3*(qd_B3F2*q_B3E1 + q_B3F2*qd_B3E1) ...
        + s13_B3*(qd_B3F1*q_B3E1 + q_B3F1*qd_B3E1));
HvS3_2 = W1_B3*qd_B3F1 + W2_B3*qd_B3F2 + W3_B3*qd_B3E1;
HvS3 = coprod(HvS3_1,J3(1,:)) + coprod(HvS3_2,J3(2,:)) + coprod(HvS3_3,J3(3,:));

EvS3   =  repmat(EvQ,1,1,nb) + HvS3 + cross(repmat(EwL,1,1,nb),rQS3,2);

%end of Velocity Vectors

%% Aerodynamic Loads
% Wind Interpolation 
ZBlNode_Y = zeros(nb,3);
ZBlNode_Z = zeros(nb,3);

ZBlNode_Y(:,1) = reshape(dot(repmat(rZO,1,1,nb) + repmat(rOQ,1,1,nb) + rQS1, repmat(Z(2,:),1,1,nb),2),nb,1);
ZBlNode_Y(:,2) = reshape(dot(repmat(rZO,1,1,nb) + repmat(rOQ,1,1,nb) + rQS2, repmat(Z(2,:),1,1,nb),2),nb,1);
ZBlNode_Y(:,3) = reshape(dot(repmat(rZO,1,1,nb) + repmat(rOQ,1,1,nb) + rQS3, repmat(Z(2,:),1,1,nb),2),nb,1);

ZBlNode_Z(:,1) = reshape(dot(repmat(rZO,1,1,nb) + repmat(rOQ,1,1,nb) + rQS1, repmat(Z(3,:),1,1,nb),2),nb,1);
ZBlNode_Z(:,2) = reshape(dot(repmat(rZO,1,1,nb) + repmat(rOQ,1,1,nb) + rQS2, repmat(Z(3,:),1,1,nb),2),nb,1);
ZBlNode_Z(:,3) = reshape(dot(repmat(rZO,1,1,nb) + repmat(rOQ,1,1,nb) + rQS3, repmat(Z(3,:),1,1,nb),2),nb,1);

% Check whether interpolation points within grid.
CheckInterpPoints(WindNom.y, WindNom.z, ZBlNode_Y, ZBlNode_Z);

Vx = zeros(nb, 3);
Vy = zeros(nb, 3);
Vz = zeros(nb, 3);

for iBd = 1:3
    for iNd = 1:nb
    Vx(iNd, iBd) = interpn([1 2 3], WindNom.y, WindNom.z, WindNom.velocity, 1, ZBlNode_Z(iNd,iBd), ZBlNode_Y(iNd,iBd));
    Vy(iNd, iBd) = interpn([1 2 3], WindNom.y, WindNom.z, WindNom.velocity, 2, ZBlNode_Z(iNd,iBd), ZBlNode_Y(iNd,iBd));
    Vz(iNd, iBd) = interpn([1 2 3], WindNom.y, WindNom.z, WindNom.velocity, 3, ZBlNode_Z(iNd,iBd), ZBlNode_Y(iNd,iBd));      
    end
end

VzB1 = coprod(reshape(Vx(:, 1),1,1,nb),Z(1,:)) + coprod(reshape(Vy(:, 1),1,1,nb),Z(3,:)) + coprod(reshape(Vz(:, 1),1,1,nb),Z(2,:));
VzB2 = coprod(reshape(Vx(:, 2),1,1,nb),Z(1,:)) + coprod(reshape(Vy(:, 2),1,1,nb),Z(3,:)) + coprod(reshape(Vz(:, 2),1,1,nb),Z(2,:));
VzB3 = coprod(reshape(Vx(:, 3),1,1,nb),Z(1,:)) + coprod(reshape(Vy(:, 3),1,1,nb),Z(3,:)) + coprod(reshape(Vz(:, 3),1,1,nb),Z(2,:));

% Blade velocities at blade local co-ordinates
UMB = zeros(nb,2,3);


UMB(:,1,1) = dot(EvS1,MB1(1,:,:),2);
UMB(:,1,2) = dot(EvS2,MB2(1,:,:),2);
UMB(:,1,3) = dot(EvS3,MB3(1,:,:),2);

if WindNom.AeroElastic
    UMB(:,2,1) = dot(EvS1,-MB1(2,:,:),2);
    UMB(:,2,2) = dot(EvS2,-MB2(2,:,:),2);       
    UMB(:,2,3) = dot(EvS3,-MB3(2,:,:),2); 
else
    UMB(:,2,1) = dot(cross(repmat(EwL,1,1,nb),rQS1,2),-MB1(2,:,:),2);
    UMB(:,2,2) = dot(cross(repmat(EwL,1,1,nb),rQS2,2),-MB2(2,:,:),2);
    UMB(:,2,3) = dot(cross(repmat(EwL,1,1,nb),rQS3,2),-MB3(2,:,:),2);
end

% Average wind vector from the field
    
AvgWndVect = mean(mean(WindNom.velocity,3),2);

if WindNom.PittandPeters
   chi0 = atan2(norm(cross(C(1,:), AvgWndVect)), dot(C(1,:), AvgWndVect)); 
else
   chi0 = 0;
end
% chi0 = 0;
% rows -> blade nodes; columns -> x, z, y directions respectively; 3rd dimension -> blade number  

% Wind speeds at blade local coordinates
VMB = zeros(nb,2,3);

VMB(:,1,1)   = dot(VzB1,MB1(1,:,:),2);          % wind velocity at blade nodes
VMB(:,2,1)   = dot(VzB1,MB1(2,:,:),2);
% VMB(:,3,1)   = dot(VzB1,MB1(3,:,:),2);        % we don't really need the z co-ordinate

VMB(:,1,2)   = dot(VzB2,MB2(1,:,:),2);
VMB(:,2,2)   = dot(VzB2,MB2(2,:,:),2);
% VMB(:,3,2)   = dot(VzB2,MB2(3,:,:),2);

VMB(:,1,3)   = dot(VzB3,MB3(1,:,:),2);
VMB(:,2,3)   = dot(VzB3,MB3(2,:,:),2);
% VMB(:,3,3)   = dot(VzB3,MB3(3,:,:),2); 
% disp(UMB)
% pause
% Load vectors

FlexBlSpn = zeros(nb,3);

FlexBlSpn(:,1) = dot(rQS1,repmat(J1(3,:),1,1,nb),2);
FlexBlSpn(:,2) = dot(rQS2,repmat(J2(3,:),1,1,nb),2);
FlexBlSpn(:,3) = dot(rQS3,repmat(J3(3,:),1,1,nb),2);


%     [px, py, Mz] = BEMT(t, q_GeAz, VMB, UMB, nb, BldSec+HubRad, FlexBlSpn, TipRad, HubRad, Cord, Twist,...
%                            BlPitch, chi0, FoilNo, Airfoils, Bld.AeroCentJ1, Bld.AeroCentJ2);
phi = reshape(Controls(5:3*nb+4), [nb 3]);

[px, py, Mz, phi] = BEMTMex(q_GeAz, VMB, UMB, phi, nb, BldSec+HubRad, FlexBlSpn, TipRad, HubRad, Cord, Bld.Twist,...
                   BlPitch, chi0, FoilNo, Airfoils, Bld.AeroCentJ1, Bld.AeroCentJ2);       

% [px, py, Mz, phi] = BEMTMex_mex(q_GeAz, VMB, UMB, phi, nb, BldSec+HubRad, FlexBlSpn, TipRad, HubRad, Cord, Bld.Twist,...
%                    BlPitch, chi0, FoilNo, Airfoils, Bld.AeroCentJ1, Bld.AeroCentJ2);       

F_AeroB1 = coprod(reshape(px(:,1),1,1,nb),MB1(1,:,:)) + coprod(reshape(py(:,1),1,1,nb),MB1(2,:,:)); 
F_AeroB2 = coprod(reshape(px(:,2),1,1,nb),MB2(1,:,:)) + coprod(reshape(py(:,2),1,1,nb),MB2(2,:,:));
F_AeroB3 = coprod(reshape(px(:,3),1,1,nb),MB3(1,:,:)) + coprod(reshape(py(:,3),1,1,nb),MB3(2,:,:));

% M_AeroB1 = coprod(reshape(Mz(:,1),1,1,nb),MB1(3,:,:));
% M_AeroB2 = coprod(reshape(Mz(:,2),1,1,nb),MB2(3,:,:));
% M_AeroB3 = coprod(reshape(Mz(:,3),1,1,nb),MB3(3,:,:));

M_AeroB1 = zeros(1,3,nb); M_AeroB2 = M_AeroB1; M_AeroB3 = M_AeroB1;

F_AeroTA1 = zeros(nt,1); F_AeroTA2 = F_AeroTA1; F_AeroTA3 = F_AeroTA1;
F_AeroT = coprod(reshape(F_AeroTA1,1,1,nt),A(1,:)) + coprod(reshape(F_AeroTA2,1,1,nt),A(2,:)) + coprod(reshape(F_AeroTA3,1,1,nt),A(3,:));

%% Hydrodynamic Loads
X  = [q_Sg, q_Sw, q_Hv, q_R, q_P, q_Y];
XD = [qd_Sg, qd_Sw, qd_Hv, qd_R, qd_P, qd_Y];


% Hydrostatic forces
C_HydroStat = [0 0       0       0            0  0; 
               0 0       0       0            0  0;
               0 0  332941       0            0  0;
               0 0       0 -4999180000         0  0;
               0 0       0       0   -4999180000  0;
               0 0       0       0            0  0];             
ff_HydroStat = [0 0 80708100 0 0 0]' - C_HydroStat*X';

% Additional damping
AddDamp = [100000 0       0   0   0     0;
            0  100000     0   0   0     0;
            0     0   130000  0   0     0;
            0     0       0   0   0     0;
            0     0       0   0   0     0;
            0     0       0   0   0 13000000];
f_AddDamp = -AddDamp*XD';

% Hyrdro Loads
F_Hydro = (ff_HydroStat(1) + f_Morison(1) + mooring_load(1)).*Z(1,:) + (ff_HydroStat(2) + f_Morison(2) + mooring_load(2)).*(-Z(3,:))...
        + (ff_HydroStat(3) + f_Morison(3) + mooring_load(3)).*Z(2,:);
M_Hydro = (ff_HydroStat(4) + f_Morison(4) + mooring_load(4)).*Z(1,:) + (ff_HydroStat(5) + f_Morison(5) + mooring_load(5)).*(-Z(3,:))...
        + (ff_HydroStat(6) + f_Morison(6) + mooring_load(6)).*Z(2,:);

%% Partial Angular Velocities

% #1. Platform

EwX_R = Z(1,:);
EwX_P = -Z(3,:);
EwX_Y = Z(2,:);

% #2. Tower Element

% All these quantities are not used
% EwF_R = Z(1,:);
% EwF_P = -Z(3,:);
% EwF_Y = Z(2,:);
% EwF_TFA1 = -do1_TFA*A(3,:);      
% EwF_TSS1 =  do1_TSS*A(1,:);   
% EwF_TFA2 = -do2_TFA*A(3,:);
% EwF_TSS2 =  do2_TSS*A(1,:);

% #3. Tower top/Base plate 

EwB_R = EwX_R;
EwB_P = EwX_P;
EwB_Y = EwX_Y;
EwB_TFA1 = -dO1_TFA*A(3,:);
EwB_TSS1 =  dO1_TSS*A(1,:);
EwB_TFA2 = -dO2_TFA*A(3,:);
EwB_TSS2 =  dO2_TSS*A(1,:);

% #4. Nacelle

EwN_R = EwB_R;
EwN_P = EwB_P;
EwN_Y = EwB_Y;
EwN_TFA1 = -dO1_TFA*A(3,:);
EwN_TSS1 =  dO1_TSS*A(1,:);
EwN_TFA2 = -dO2_TFA*A(3,:);
EwN_TSS2 =  dO2_TSS*A(1,:);
EwN_yaw  =  D(2,:);

% #5. Low speed shaft

EwL_R = EwN_R;
EwL_P = EwN_P;
EwL_Y = EwN_Y;
EwL_TFA1 = -dO1_TFA*A(3,:);
EwL_TSS1 =  dO1_TSS*A(1,:);
EwL_TFA2 = -dO2_TFA*A(3,:);
EwL_TSS2 =  dO2_TSS*A(1,:);
EwL_yaw  =  D(2,:);
EwL_GeAz  =  C(1,:);
EwL_DrTr  =  C(1,:);

% #6. Blade 1 Element

% EwM1_R = repmat(EwN_R,1,1,nb);
% EwM1_P = repmat(EwN_P,1,1,nb);
% EwM1_Y = repmat(EwN_Y,1,1,nb);
% EwM1_TFA1 = repmat(EwL_TFA1,1,1,nb);
% EwM1_TSS1 = repmat(EwL_TSS1,1,1,nb);
% EwM1_TFA2 = repmat(EwL_TFA2,1,1,nb);
% EwM1_TSS2 = repmat(EwL_TSS2,1,1,nb);
% EwM1_yaw  =  repmat(EwL_yaw,1,1,nb);
% EwM1_GeAz  = repmat(EwL_GeAz,1,1,nb);
% EwM1_DrTr  = repmat(EwL_DrTr,1,1,nb);
EwM1_B1F1  =  coprod(-dW1_B1,J1(1,:)) + coprod(dO1_B1,J1(2,:));
EwM1_B1E1  =  coprod(-dW3_B1,J1(1,:)) + coprod(dO3_B1,J1(2,:));
EwM1_B1F2  =  coprod(-dW2_B1,J1(1,:)) + coprod(dO2_B1,J1(2,:));

% #7. Blade 2 Element

% EwM2_R = repmat(EwN_R,1,1,nb);
% EwM2_P = repmat(EwN_P,1,1,nb);
% EwM2_Y = repmat(EwN_Y,1,1,nb);
% EwM2_TFA1 = repmat(EwL_TFA1,1,1,nb);
% EwM2_TSS1 = repmat(EwL_TSS1,1,1,nb);
% EwM2_TFA2 = repmat(EwL_TFA2,1,1,nb);
% EwM2_TSS2 = repmat(EwL_TSS2,1,1,nb);
% EwM2_yaw  =  repmat(EwL_yaw,1,1,nb);
% EwM2_GeAz  = repmat(EwL_GeAz,1,1,nb);
% EwM2_DrTr  = repmat(EwL_DrTr,1,1,nb);
EwM2_B2F1  =  coprod(-dW1_B2,J2(1,:)) + coprod(dO1_B2,J2(2,:));
EwM2_B2E1  =  coprod(-dW3_B2,J2(1,:)) + coprod(dO3_B2,J2(2,:));
EwM2_B2F2  =  coprod(-dW2_B2,J2(1,:)) + coprod(dO2_B2,J2(2,:));

% #8. Blade 3 Element

% EwM3_R = repmat(EwN_R,1,1,nb);
% EwM3_P = repmat(EwN_P,1,1,nb);
% EwM3_Y = repmat(EwN_Y,1,1,nb);
% EwM3_TFA1 = repmat(EwL_TFA1,1,1,nb);
% EwM3_TSS1 = repmat(EwL_TSS1,1,1,nb);
% EwM3_TFA2 = repmat(EwL_TFA2,1,1,nb);
% EwM3_TSS2 =  repmat(EwL_TSS2,1,1,nb);
% EwM3_yaw  =  repmat(EwL_yaw,1,1,nb);
% EwM3_GeAz  = repmat(EwL_GeAz,1,1,nb);
% EwM3_DrTr  = repmat(EwL_DrTr,1,1,nb);
EwM3_B3F1  =  coprod(-dW1_B3,J3(1,:)) + coprod(dO1_B3,J3(2,:));
EwM3_B3E1  =  coprod(-dW3_B3,J3(1,:)) + coprod(dO3_B3,J3(2,:));
EwM3_B3F2  =  coprod(-dW2_B3,J3(1,:)) + coprod(dO2_B3,J3(2,:));


% #9. High speed shaft/Generator

EwG_R = Z(1,:);
EwG_P = -Z(3,:);
EwG_Y = Z(2,:);
EwG_TFA1 = -dO1_TFA*A(3,:);
EwG_TSS1 =  dO1_TSS*A(1,:);
EwG_TFA2 = -dO2_TFA*A(3,:);
EwG_TSS2 =  dO2_TSS*A(1,:);
EwG_yaw  =  D(2,:);
EwG_GeAz = GenDir*GBRatio*C(1,:);

% end of Partial Angular Velocities

%% Partial Linear Velocities

% #1. Platform

EvZ_Sg = Z(1,:);
EvZ_Sw = -Z(3,:);
EvZ_Hv = Z(2,:);

% #2. Platform C.M
EvY_Sg = EvZ_Sg;
EvY_Sw = EvZ_Sw;
EvY_Hv = EvZ_Hv;
EvY_R  = cross(EwX_R,rZY);
EvY_P  = cross(EwX_P,rZY);
EvY_Y  = cross(EwX_Y,rZY);

% #3. Tower

EvF_Sg = repmat(Z(1,:),[1 1 nt]);
EvF_Sw = repmat(-Z(3,:),[1 1 nt]);
EvF_Hv = repmat(Z(2,:),[1 1 nt]);
EvF_R  = cross(repmat(EwX_R,[1 1 nt]),rZT,2);
EvF_P  = cross(repmat(EwX_P,[1 1 nt]),rZT,2);
EvF_Y  = cross(repmat(EwX_Y,[1 1 nt]),rZT,2);
EvF_TFA1 = coprod(O1_TFA,A(1,:)) - coprod((s11_TFA*q_TFA1 + s12_TFA*q_TFA2),A(2,:)); 
EvF_TSS1 = coprod(O1_TSS,A(3,:)) - coprod((s11_TSS*q_TSS1 + s12_TSS*q_TSS2),A(2,:));
EvF_TFA2 = coprod(O2_TFA,A(1,:)) - coprod((s22_TFA*q_TFA2 + s12_TFA*q_TFA1),A(2,:));
EvF_TSS2 = coprod(O2_TSS,A(3,:)) - coprod((s22_TSS*q_TSS2 + s12_TSS*q_TSS1),A(2,:));

% #4. Tower top/Base plate

EvO_Sg = EvY_Sg;
EvO_Sw = EvY_Sw;
EvO_Hv = EvY_Hv;
EvO_R  = cross(EwX_R,rZO);
EvO_P  = cross(EwX_P,rZO);
EvO_Y  = cross(EwX_Y,rZO);
EvO_TFA1 = A(1,:) - (S11_TFA*q_TFA1 + S12_TFA*q_TFA2)*A(2,:);
EvO_TSS1 = A(3,:) - (S11_TSS*q_TSS1 + S12_TSS*q_TSS2)*A(2,:);
EvO_TFA2 = A(1,:) - (S22_TFA*q_TFA2 + S12_TFA*q_TFA1)*A(2,:);
EvO_TSS2 = A(3,:) - (S22_TSS*q_TSS2 + S12_TSS*q_TSS1)*A(2,:);

% #5. Nacelle
EvU_Sg = EvO_Sg;
EvU_Sw = EvO_Sw;
EvU_Hv = EvO_Hv;
EvU_R  = EvO_R + cross(EwN_R,rOU);
EvU_P  = EvO_P + cross(EwN_P,rOU);
EvU_Y  = EvO_Y + cross(EwN_Y,rOU);
EvU_TFA1  = EvO_TFA1 + cross(EwN_TFA1,rOU);
EvU_TSS1  = EvO_TSS1 + cross(EwN_TSS1,rOU);
EvU_TFA2  = EvO_TFA2 + cross(EwN_TFA2,rOU);
EvU_TSS2  = EvO_TSS2 + cross(EwN_TSS2,rOU);
EvU_yaw   = cross(EwN_yaw,rOU);

% #6. Appex of conning angle

EvQ_Sg = EvO_Sg;
EvQ_Sw = EvO_Sw;
EvQ_Hv = EvO_Hv;
EvQ_R  = EvO_R + cross(EwN_R,rOQ);
EvQ_P  = EvO_P + cross(EwN_P,rOQ);
EvQ_Y  = EvO_Y + cross(EwN_Y,rOQ);
EvQ_TFA1  = EvO_TFA1 + cross(EwN_TFA1,rOQ);
EvQ_TSS1  = EvO_TSS1 + cross(EwN_TSS1,rOQ);
EvQ_TFA2  = EvO_TFA2 + cross(EwN_TFA2,rOQ);
EvQ_TSS2  = EvO_TSS2 + cross(EwN_TSS2,rOQ);
EvQ_yaw   = cross(EwN_yaw,rOQ);
% EvQ_GeAz  = cross(EwN_GeAz,rOQ);
% EvQ_DrTr  = cross(EwN_DrTr,rOQ);

% #7. Hub C.M

EvC_Sg = EvQ_Sg;
EvC_Sw = EvQ_Sw;
EvC_Hv = EvQ_Hv;
EvC_R  = EvQ_R + cross(EwL_R,rQC);
EvC_P  = EvQ_P + cross(EwL_P,rQC);
EvC_Y  = EvQ_Y + cross(EwL_Y,rQC);
EvC_TFA1  = EvQ_TFA1 + cross(EwL_TFA1,rQC);
EvC_TSS1  = EvQ_TSS1 + cross(EwL_TSS1,rQC);
EvC_TFA2  = EvQ_TFA2 + cross(EwL_TFA2,rQC);
EvC_TSS2  = EvQ_TSS2 + cross(EwL_TSS2,rQC);
EvC_yaw   = EvQ_yaw  + cross(EwL_yaw,rQC);
EvC_GeAz  = cross(EwL_GeAz,rQC);
EvC_DrTr  = cross(EwL_DrTr,rQC);

% #8. Blade 1

EvS1_Sg = repmat(EvQ_Sg,1,1,nb);
EvS1_Sw = repmat(EvQ_Sw,1,1,nb);
EvS1_Hv = repmat(EvQ_Hv,1,1,nb);
EvS1_R  = repmat(EvQ_R,1,1,nb) + cross(repmat(EwL_R,1,1,nb),rQS1,2);
EvS1_P  = repmat(EvQ_P,1,1,nb) + cross(repmat(EwL_P,1,1,nb),rQS1,2);
EvS1_Y  = repmat(EvQ_Y,1,1,nb) + cross(repmat(EwL_Y,1,1,nb),rQS1,2);
EvS1_TFA1  = repmat(EvQ_TFA1,1,1,nb) + cross(repmat(EwL_TFA1,1,1,nb),rQS1,2);
EvS1_TSS1  = repmat(EvQ_TSS1,1,1,nb) + cross(repmat(EwL_TSS1,1,1,nb),rQS1,2);
EvS1_TFA2  = repmat(EvQ_TFA2,1,1,nb) + cross(repmat(EwL_TFA2,1,1,nb),rQS1,2);
EvS1_TSS2  = repmat(EvQ_TSS2,1,1,nb) + cross(repmat(EwL_TSS2,1,1,nb),rQS1,2);
EvS1_yaw   = repmat(EvQ_yaw,1,1,nb)  + cross(repmat(EwL_yaw,1,1,nb),rQS1,2);
EvS1_GeAz  = cross(repmat(EwL_GeAz,1,1,nb),rQS1,2);
EvS1_DrTr  = cross(repmat(EwL_DrTr,1,1,nb),rQS1,2);
EvS1_B1F1  = coprod(O1_B1,J1(1,:)) + coprod(W1_B1,J1(2,:)) - coprod((s11_B1*q_B1F1 + s12_B1*q_B1F2 + s13_B1*q_B1E1),J1(3,:)); 
EvS1_B1E1  = coprod(O3_B1,J1(1,:)) + coprod(W3_B1,J1(2,:)) - coprod((s33_B1*q_B1E1 + s23_B1*q_B1F2 + s13_B1*q_B1F1),J1(3,:));
EvS1_B1F2  = coprod(O2_B1,J1(1,:)) + coprod(W2_B1,J1(2,:)) - coprod((s22_B1*q_B1F2 + s12_B1*q_B1F1 + s23_B1*q_B1E1),J1(3,:));

% Blade 2

EvS2_Sg = repmat(EvQ_Sg,1,1,nb);
EvS2_Sw = repmat(EvQ_Sw,1,1,nb);
EvS2_Hv = repmat(EvQ_Hv,1,1,nb);
EvS2_R  = repmat(EvQ_R,1,1,nb) + cross(repmat(EwL_R,1,1,nb),rQS2,2);
EvS2_P  = repmat(EvQ_P,1,1,nb) + cross(repmat(EwL_P,1,1,nb),rQS2,2);
EvS2_Y  = repmat(EvQ_Y,1,1,nb) + cross(repmat(EwL_Y,1,1,nb),rQS2,2);
EvS2_TFA1  = repmat(EvQ_TFA1,1,1,nb) + cross(repmat(EwL_TFA1,1,1,nb),rQS2,2);
EvS2_TSS1  = repmat(EvQ_TSS1,1,1,nb) + cross(repmat(EwL_TSS1,1,1,nb),rQS2,2);
EvS2_TFA2  = repmat(EvQ_TFA2,1,1,nb) + cross(repmat(EwL_TFA2,1,1,nb),rQS2,2);
EvS2_TSS2  = repmat(EvQ_TSS2,1,1,nb) + cross(repmat(EwL_TSS2,1,1,nb),rQS2,2);
EvS2_yaw   = repmat(EvQ_yaw,1,1,nb)  + cross(repmat(EwL_yaw,1,1,nb),rQS2,2);
EvS2_GeAz  = cross(repmat(EwL_GeAz,1,1,nb),rQS2,2);
EvS2_DrTr  = cross(repmat(EwL_DrTr,1,1,nb),rQS2,2);
EvS2_B2F1  = coprod(O1_B2,J2(1,:)) + coprod(W1_B2,J2(2,:)) - coprod((s11_B2*q_B2F1 + s12_B2*q_B2F2 + s13_B2*q_B2E1),J2(3,:));
EvS2_B2E1  = coprod(O3_B2,J2(1,:)) + coprod(W3_B2,J2(2,:)) - coprod((s33_B2*q_B2E1 + s23_B2*q_B2F2 + s13_B2*q_B2F1),J2(3,:));
EvS2_B2F2  = coprod(O2_B2,J2(1,:)) + coprod(W2_B2,J2(2,:)) - coprod((s22_B2*q_B2F2 + s12_B2*q_B2F1 + s23_B2*q_B2E1),J2(3,:));

% Blade 3

EvS3_Sg = repmat(EvQ_Sg,1,1,nb);
EvS3_Sw = repmat(EvQ_Sw,1,1,nb);
EvS3_Hv = repmat(EvQ_Hv,1,1,nb);
EvS3_R  = repmat(EvQ_R,1,1,nb) + cross(repmat(EwL_R,1,1,nb),rQS3,2);
EvS3_P  = repmat(EvQ_P,1,1,nb) + cross(repmat(EwL_P,1,1,nb),rQS3,2);
EvS3_Y  = repmat(EvQ_Y,1,1,nb) + cross(repmat(EwL_Y,1,1,nb),rQS3,2);
EvS3_TFA1  = repmat(EvQ_TFA1,1,1,nb) + cross(repmat(EwL_TFA1,1,1,nb),rQS3,2);
EvS3_TSS1  = repmat(EvQ_TSS1,1,1,nb) + cross(repmat(EwL_TSS1,1,1,nb),rQS3,2);
EvS3_TFA2  = repmat(EvQ_TFA2,1,1,nb) + cross(repmat(EwL_TFA2,1,1,nb),rQS3,2);
EvS3_TSS2  = repmat(EvQ_TSS2,1,1,nb) + cross(repmat(EwL_TSS2,1,1,nb),rQS3,2);
EvS3_yaw   = repmat(EvQ_yaw,1,1,nb)  + cross(repmat(EwL_yaw,1,1,nb),rQS3,2);
EvS3_GeAz  = cross(repmat(EwL_GeAz,1,1,nb),rQS3,2);
EvS3_DrTr  = cross(repmat(EwL_DrTr,1,1,nb),rQS3,2);
EvS3_B3F1  = coprod(O1_B3,J3(1,:)) + coprod(W1_B3,J3(2,:)) - coprod((s11_B3*q_B3F1 + s12_B3*q_B3F2 + s13_B3*q_B3E1),J3(3,:));
EvS3_B3E1  = coprod(O3_B3,J3(1,:)) + coprod(W3_B3,J3(2,:)) - coprod((s33_B3*q_B3E1 + s23_B3*q_B3F2 + s13_B3*q_B3F1),J3(3,:));
EvS3_B3F2  = coprod(O2_B3,J3(1,:)) + coprod(W2_B3,J3(2,:)) - coprod((s22_B3*q_B3F2 + s12_B3*q_B3F1 + s23_B3*q_B3E1),J3(3,:));


% end of Partial Linear Velocities

%% Angular Accelerations

% #1. Platform
% all zeros

% #2. Tower

% I think not required

% #3. Tower top/Base plate

dEwB_TFA1 = cross(EwX,EwB_TFA1);
dEwB_TSS1 = cross(EwX,EwB_TSS1);
dEwB_TFA2 = cross(EwX,EwB_TFA2);
dEwB_TSS2 = cross(EwX,EwB_TSS2);

% #4. Nacelle

dEwN_TFA1 = dEwB_TFA1;
dEwN_TSS1 = dEwB_TSS1;
dEwN_TFA2 = dEwB_TFA2;
dEwN_TSS2 = dEwB_TSS2;
dEwN_yaw  = cross(EwB,EwN_yaw);

dEwN = dEwN_TFA1*qd_TFA1 + dEwN_TSS1*qd_TSS1 + dEwN_TFA2*qd_TFA2 + dEwN_TSS2*qd_TSS2 + dEwN_yaw*qd_yaw;

% #5. Low Speed Shaft

dEwL_TFA1 = dEwN_TFA1;
dEwL_TSS1 = dEwN_TSS1;
dEwL_TFA2 = dEwN_TFA2;
dEwL_TSS2 = dEwN_TSS2;
dEwL_yaw  = dEwN_yaw;
dEwL_GeAz = cross(EwN,EwL_GeAz);
dEwL_DrTr = cross(EwN,EwL_DrTr);

dEwL = dEwL_TFA1*qd_TFA1 + dEwL_TSS1*qd_TSS1 + dEwL_TFA2*qd_TFA2 + dEwL_TSS2*qd_TSS2 + dEwL_yaw*qd_yaw + dEwL_GeAz*qd_GeAz + dEwL_DrTr*qd_DrTr;

% #6. Blade 1

% dEwM1_TFA1 = dEwL_TFA1;
% dEwM1_TSS1 = dEwL_TSS1;
% dEwM1_TFA2 = dEwL_TFA2;
% dEwM1_TSS2 = dEwL_TSS2;
% dEwM1_yaw  = dEwL_yaw;
% dEwM1_GeAz = dEwL_GeAz;
% dEwM1_DrTr = dEwL_DrTr;
% dEwM1_B1F1 = cross(EwL,EwM1_B1F1);
% dEwM1_B1E1 = cross(EwL,EwM1_B1E1);
% dEwM1_B1F2 = cross(EwL,EwM1_B1F2);

% #7. High Speed Shaft/Generator

dEwG_TFA1 = dEwN_TFA1;
dEwG_TSS1 = dEwN_TSS1;
dEwG_TFA2 = dEwN_TFA2;
dEwG_TSS2 = dEwN_TSS2;
dEwG_yaw  = dEwN_yaw;
dEwG_GeAz = cross(EwN,EwG_GeAz);

dEwG      = dEwG_TFA1*qd_TFA1 + dEwG_TSS1*qd_TSS1 + dEwG_TFA2*qd_TFA2 + dEwG_TSS2*qd_TSS2 + dEwG_yaw*qd_yaw + dEwG_GeAz*qd_GeAz;
% end of Angular Accelerations

%% Linear accelerations

% #2. Platform C.M

dEvY_R = cross(EwX_R,cross(EwX,rZY)); 
dEvY_P = cross(EwX_P,cross(EwX,rZY));
dEvY_Y = cross(EwX_Y,cross(EwX,rZY));

dEvY = dEvY_R*qd_R + dEvY_P*qd_P + dEvY_Y*qd_Y;

% #3. Tower

dEvF_R = cross(repmat(EwX_R,[1, 1, nt]), XvT + cross(repmat(EwX,[1, 1, nt]),rZT,2), 2);
dEvF_P = cross(repmat(EwX_P,[1, 1, nt]), XvT + cross(repmat(EwX,[1, 1, nt]),rZT,2), 2);
dEvF_Y = cross(repmat(EwX_Y,[1, 1, nt]), XvT + cross(repmat(EwX,[1, 1, nt]),rZT,2), 2);
dEvF_TFA1 = -coprod((s11_TFA*qd_TFA1+s12_TFA*qd_TFA2),A(2,:)) + cross(repmat(EwX,[1, 1, nt]), EvF_TFA1,2);     
dEvF_TSS1 = -coprod((s11_TSS*qd_TSS1+s12_TSS*qd_TSS2),A(2,:)) + cross(repmat(EwX,[1, 1, nt]), EvF_TSS1,2);
dEvF_TFA2 = -coprod((s22_TFA*qd_TFA2+s12_TFA*qd_TFA1),A(2,:)) + cross(repmat(EwX,[1, 1, nt]), EvF_TFA2,2);
dEvF_TSS2 = -coprod((s22_TSS*qd_TSS2+s12_TSS*qd_TSS1),A(2,:)) + cross(repmat(EwX,[1, 1, nt]), EvF_TSS2,2);

dEvF = dEvF_R*qd_R + dEvF_P*qd_P + dEvF_Y*qd_Y + dEvF_TFA1*qd_TFA1 + dEvF_TSS1*qd_TSS1 + dEvF_TFA2*qd_TFA2 + dEvF_TSS2*qd_TSS2;

% #4. Tower top/Base plate

dEvO_R = cross(EwX_R, XvO + cross(EwX,rZO));
dEvO_P = cross(EwX_P, XvO + cross(EwX,rZO));
dEvO_Y = cross(EwX_Y, XvO + cross(EwX,rZO));
dEvO_TFA1 = -(S11_TFA*qd_TFA1 + S12_TFA*qd_TFA2)*A(2,:) + cross(EwX, EvO_TFA1);
dEvO_TSS1 = -(S11_TSS*qd_TSS1 + S12_TSS*qd_TSS2)*A(2,:) + cross(EwX, EvO_TSS1);
dEvO_TFA2 = -(S22_TFA*qd_TFA2 + S12_TFA*qd_TFA1)*A(2,:) + cross(EwX, EvO_TFA2);
dEvO_TSS2 = -(S22_TSS*qd_TSS2 + S12_TSS*qd_TSS1)*A(2,:) + cross(EwX, EvO_TSS2);

dEvO = dEvO_R*qd_R + dEvO_P*qd_P + dEvO_Y*qd_Y + dEvO_TFA1*qd_TFA1 + dEvO_TSS1*qd_TSS1 + dEvO_TFA2*qd_TFA2 + dEvO_TSS2*qd_TSS2;

% #5. Nacelle

dEvU_R = dEvO_R + cross(EwN_R ,cross(EwN,rOU));
dEvU_P = dEvO_P + cross(EwN_P ,cross(EwN,rOU));
dEvU_Y = dEvO_Y + cross(EwN_Y ,cross(EwN,rOU));
dEvU_TFA1 = dEvO_TFA1 + cross(dEwN_TFA1,rOU) + cross(EwN_TFA1,cross(EwN,rOU));
dEvU_TSS1 = dEvO_TSS1 + cross(dEwN_TSS1,rOU) + cross(EwN_TSS1,cross(EwN,rOU));
dEvU_TFA2 = dEvO_TFA2 + cross(dEwN_TFA2,rOU) + cross(EwN_TFA2,cross(EwN,rOU));
dEvU_TSS2 = dEvO_TSS2 + cross(dEwN_TSS2,rOU) + cross(EwN_TSS2,cross(EwN,rOU));
dEvU_yaw  = cross(dEwN_yaw,rOU) + cross(EwN_yaw,cross(EwN,rOU));

dEvU = dEvU_R*qd_R + dEvU_P*qd_P + dEvU_Y*qd_Y + dEvU_TFA1*qd_TFA1 + dEvU_TSS1*qd_TSS1 + dEvU_TFA2*qd_TFA2 + dEvU_TSS2*qd_TSS2 + dEvU_yaw*qd_yaw;

% #6. Appex of Conning angle

dEvQ_R = dEvO_R + cross(EwN_R ,cross(EwN,rOQ));
dEvQ_P = dEvO_P + cross(EwN_P ,cross(EwN,rOQ));
dEvQ_Y = dEvO_Y + cross(EwN_Y ,cross(EwN,rOQ));
dEvQ_TFA1 = dEvO_TFA1 + cross(dEwN_TFA1,rOQ) + cross(EwN_TFA1,cross(EwN,rOQ));
dEvQ_TSS1 = dEvO_TSS1 + cross(dEwN_TSS1,rOQ) + cross(EwN_TSS1,cross(EwN,rOQ));
dEvQ_TFA2 = dEvO_TFA2 + cross(dEwN_TFA2,rOQ) + cross(EwN_TFA2,cross(EwN,rOQ));
dEvQ_TSS2 = dEvO_TSS2 + cross(dEwN_TSS2,rOQ) + cross(EwN_TSS2,cross(EwN,rOQ));
dEvQ_yaw  = cross(dEwN_yaw,rOQ) + cross(EwN_yaw,cross(EwN,rOQ));
% dEvQ_GeAz = cross(dEwL_GeAz,rOQ) + cross(EwL_GeAz,cross(EwL,rOQ));
% dEvQ_DrTr = cross(dEwL_DrTr,rOQ) + cross(EwL_DrTr,cross(EwL,rOQ));

% #7. Hub C.M

dEvC_R = dEvQ_R + cross(EwL_R ,cross(EwL,rQC));
dEvC_P = dEvQ_P + cross(EwL_P ,cross(EwL,rQC));
dEvC_Y = dEvQ_Y + cross(EwL_Y ,cross(EwL,rQC));
dEvC_TFA1 = dEvQ_TFA1 + cross(dEwL_TFA1,rQC) + cross(EwL_TFA1,cross(EwL,rQC));
dEvC_TSS1 = dEvQ_TSS1 + cross(dEwL_TSS1,rQC) + cross(EwL_TSS1,cross(EwL,rQC));
dEvC_TFA2 = dEvQ_TFA2 + cross(dEwL_TFA2,rQC) + cross(EwL_TFA2,cross(EwL,rQC));
dEvC_TSS2 = dEvQ_TSS2 + cross(dEwL_TSS2,rQC) + cross(EwL_TSS2,cross(EwL,rQC));
dEvC_yaw  = dEvQ_yaw  + cross(dEwL_yaw,rQC)  + cross(EwL_yaw,cross(EwL,rQC));
dEvC_GeAz = cross(dEwL_GeAz,rQC) + cross(EwL_GeAz,cross(EwL,rQC)); 
dEvC_DrTr = cross(dEwL_DrTr,rQC) + cross(EwL_DrTr,cross(EwL,rQC));

dEvC = dEvC_R*qd_R + dEvC_P*qd_P + dEvC_Y*qd_Y + dEvC_TFA1*qd_TFA1 + dEvC_TSS1*qd_TSS1 + dEvC_TFA2*qd_TFA2 + dEvC_TSS2*qd_TSS2 ...
       + dEvC_yaw*qd_yaw + dEvC_GeAz*qd_GeAz + dEvC_DrTr*qd_DrTr;

% #8. Blade 1

dEvS1_R = repmat(dEvQ_R,1,1,nb) + cross(repmat(EwL_R,[1,1,nb]), HvS1 + cross(repmat(EwL,[1,1,nb]), rQS1, 2), 2);
dEvS1_P = repmat(dEvQ_P,1,1,nb) + cross(repmat(EwL_P,[1,1,nb]), HvS1 + cross(repmat(EwL,[1,1,nb]), rQS1, 2), 2);
dEvS1_Y = repmat(dEvQ_Y,1,1,nb) + cross(repmat(EwL_Y,[1,1,nb]), HvS1 + cross(repmat(EwL,[1,1,nb]), rQS1, 2), 2);
dEvS1_TFA1 = repmat(dEvQ_TFA1,1,1,nb) + cross(repmat(dEwL_TFA1,[1,1,nb]),rQS1,2) + cross(repmat(EwL_TFA1,[1,1,nb]), HvS1 + cross(repmat(EwL,[1,1,nb]), rQS1, 2), 2);
dEvS1_TSS1 = repmat(dEvQ_TSS1,1,1,nb) + cross(repmat(dEwL_TSS1,[1,1,nb]),rQS1,2) + cross(repmat(EwL_TSS1,[1,1,nb]), HvS1 + cross(repmat(EwL,[1,1,nb]), rQS1, 2), 2);
dEvS1_TFA2 = repmat(dEvQ_TFA2,1,1,nb) + cross(repmat(dEwL_TFA2,[1,1,nb]),rQS1,2) + cross(repmat(EwL_TFA2,[1,1,nb]), HvS1 + cross(repmat(EwL,[1,1,nb]), rQS1, 2), 2);
dEvS1_TSS2 = repmat(dEvQ_TSS2,1,1,nb) + cross(repmat(dEwL_TSS2,[1,1,nb]),rQS1,2) + cross(repmat(EwL_TSS2,[1,1,nb]), HvS1 + cross(repmat(EwL,[1,1,nb]), rQS1, 2), 2);
dEvS1_yaw  = repmat(dEvQ_yaw,1,1,nb)  + cross(repmat(dEwL_yaw,[1,1,nb]),rQS1,2)  + cross(repmat(EwL_yaw,[1,1,nb]),  HvS1 + cross(repmat(EwL,[1,1,nb]), rQS1, 2), 2);
dEvS1_GeAz = cross(repmat(dEwL_GeAz,[1,1,nb]),rQS1,2) + cross(repmat(EwL_GeAz,[1,1,nb]), HvS1 + cross(repmat(EwL,[1,1,nb]), rQS1, 2), 2);
dEvS1_DrTr = cross(repmat(dEwL_DrTr,[1,1,nb]),rQS1,2) + cross(repmat(EwL_DrTr,[1,1,nb]), HvS1 + cross(repmat(EwL,[1,1,nb]), rQS1, 2), 2);
dEvS1_B1F1 = -coprod((s11_B1*qd_B1F1 + s12_B1*qd_B1F2 + s13_B1*qd_B1E1),J1(3,:)) + cross(repmat(EwL,1,1,nb),EvS1_B1F1,2);        % no need to reshape
dEvS1_B1E1 = -coprod((s33_B1*qd_B1E1 + s23_B1*qd_B1F2 + s13_B1*qd_B1F1),J1(3,:)) + cross(repmat(EwL,1,1,nb),EvS1_B1E1,2);
dEvS1_B1F2 = -coprod((s22_B1*qd_B1F2 + s12_B1*qd_B1F1 + s23_B1*qd_B1E1),J1(3,:)) + cross(repmat(EwL,1,1,nb),EvS1_B1F2,2);

dEvS1 = dEvS1_R*qd_R + dEvS1_P*qd_P + dEvS1_Y*qd_Y + dEvS1_TFA1*qd_TFA1 + dEvS1_TSS1*qd_TSS1 + dEvS1_TFA2*qd_TFA2 + dEvS1_TSS2*qd_TSS2 ...
      + dEvS1_yaw*qd_yaw + dEvS1_GeAz*qd_GeAz + dEvS1_DrTr*qd_DrTr + dEvS1_B1F1*qd_B1F1 + dEvS1_B1E1*qd_B1E1 + dEvS1_B1F2*qd_B1F2;

% #9. Blade 2

dEvS2_R = repmat(dEvQ_R,1,1,nb) + cross(repmat(EwL_R,[1,1,nb]), HvS2 + cross(repmat(EwL,[1,1,nb]), rQS2, 2), 2);
dEvS2_P = repmat(dEvQ_P,1,1,nb) + cross(repmat(EwL_P,[1,1,nb]), HvS2 + cross(repmat(EwL,[1,1,nb]), rQS2, 2), 2);
dEvS2_Y = repmat(dEvQ_Y,1,1,nb) + cross(repmat(EwL_Y,[1,1,nb]), HvS2 + cross(repmat(EwL,[1,1,nb]), rQS2, 2), 2);
dEvS2_TFA1 = repmat(dEvQ_TFA1,1,1,nb) + cross(repmat(dEwL_TFA1,[1,1,nb]),rQS2,2) + cross(repmat(EwL_TFA1,[1,1,nb]), HvS2 + cross(repmat(EwL,[1,1,nb]), rQS2, 2), 2);
dEvS2_TSS1 = repmat(dEvQ_TSS1,1,1,nb) + cross(repmat(dEwL_TSS1,[1,1,nb]),rQS2,2) + cross(repmat(EwL_TSS1,[1,1,nb]), HvS2 + cross(repmat(EwL,[1,1,nb]), rQS2, 2), 2);
dEvS2_TFA2 = repmat(dEvQ_TFA2,1,1,nb) + cross(repmat(dEwL_TFA2,[1,1,nb]),rQS2,2) + cross(repmat(EwL_TFA2,[1,1,nb]), HvS2 + cross(repmat(EwL,[1,1,nb]), rQS2, 2), 2);
dEvS2_TSS2 = repmat(dEvQ_TSS2,1,1,nb) + cross(repmat(dEwL_TSS2,[1,1,nb]),rQS2,2) + cross(repmat(EwL_TSS2,[1,1,nb]), HvS2 + cross(repmat(EwL,[1,1,nb]), rQS2, 2), 2);
dEvS2_yaw  = repmat(dEvQ_yaw,1,1,nb)  + cross(repmat(dEwL_yaw,[1,1,nb]),rQS2,2)  + cross(repmat(EwL_yaw,[1,1,nb]),  HvS2 + cross(repmat(EwL,[1,1,nb]), rQS2, 2), 2);
dEvS2_GeAz = cross(repmat(dEwL_GeAz,[1,1,nb]),rQS2,2) + cross(repmat(EwL_GeAz,[1,1,nb]), HvS2 + cross(repmat(EwL,[1,1,nb]), rQS2, 2), 2);
dEvS2_DrTr = cross(repmat(dEwL_DrTr,[1,1,nb]),rQS2,2) + cross(repmat(EwL_DrTr,[1,1,nb]), HvS2 + cross(repmat(EwL,[1,1,nb]), rQS2, 2), 2);
dEvS2_B2F1 = -coprod((s11_B2*qd_B2F1 + s12_B2*qd_B2F2 + s13_B2*qd_B2E1),J2(3,:)) + cross(repmat(EwL,1,1,nb),EvS2_B2F1,2);
dEvS2_B2E1 = -coprod((s33_B2*qd_B2E1 + s23_B2*qd_B2F2 + s13_B2*qd_B2F1),J2(3,:)) + cross(repmat(EwL,1,1,nb),EvS2_B2E1,2);
dEvS2_B2F2 = -coprod((s22_B2*qd_B2F2 + s12_B2*qd_B2F1 + s23_B2*qd_B2E1),J2(3,:)) + cross(repmat(EwL,1,1,nb),EvS2_B2F2,2);

dEvS2 = dEvS2_R*qd_R + dEvS2_P*qd_P + dEvS2_Y*qd_Y + dEvS2_TFA1*qd_TFA1 + dEvS2_TSS1*qd_TSS1 + dEvS2_TFA2*qd_TFA2 + dEvS2_TSS2*qd_TSS2 ...
      + dEvS2_yaw*qd_yaw + dEvS2_GeAz*qd_GeAz + dEvS2_DrTr*qd_DrTr + dEvS2_B2F1*qd_B2F1 + dEvS2_B2E1*qd_B2E1 + dEvS2_B2F2*qd_B2F2;

% #10. Blade 3

dEvS3_R = repmat(dEvQ_R,1,1,nb) + cross(repmat(EwL_R,[1,1,nb]), HvS3 + cross(repmat(EwL,[1,1,nb]), rQS3, 2), 2);
dEvS3_P = repmat(dEvQ_P,1,1,nb) + cross(repmat(EwL_P,[1,1,nb]), HvS3 + cross(repmat(EwL,[1,1,nb]), rQS3, 2), 2);
dEvS3_Y = repmat(dEvQ_Y,1,1,nb) + cross(repmat(EwL_Y,[1,1,nb]), HvS3 + cross(repmat(EwL,[1,1,nb]), rQS3, 2), 2);
dEvS3_TFA1 = repmat(dEvQ_TFA1,1,1,nb) + cross(repmat(dEwL_TFA1,[1,1,nb]),rQS3,2) + cross(repmat(EwL_TFA1,[1,1,nb]), HvS3 + cross(repmat(EwL,[1,1,nb]), rQS3, 2), 2);
dEvS3_TSS1 = repmat(dEvQ_TSS1,1,1,nb) + cross(repmat(dEwL_TSS1,[1,1,nb]),rQS3,2) + cross(repmat(EwL_TSS1,[1,1,nb]), HvS3 + cross(repmat(EwL,[1,1,nb]), rQS3, 2), 2);
dEvS3_TFA2 = repmat(dEvQ_TFA2,1,1,nb) + cross(repmat(dEwL_TFA2,[1,1,nb]),rQS3,2) + cross(repmat(EwL_TFA2,[1,1,nb]), HvS3 + cross(repmat(EwL,[1,1,nb]), rQS3, 2), 2);
dEvS3_TSS2 = repmat(dEvQ_TSS2,1,1,nb) + cross(repmat(dEwL_TSS2,[1,1,nb]),rQS3,2) + cross(repmat(EwL_TSS2,[1,1,nb]), HvS3 + cross(repmat(EwL,[1,1,nb]), rQS3, 2), 2);
dEvS3_yaw  = repmat(dEvQ_yaw,1,1,nb)  + cross(repmat(dEwL_yaw,[1,1,nb]),rQS3,2)  + cross(repmat(EwL_yaw,[1,1,nb]),  HvS3 + cross(repmat(EwL,[1,1,nb]), rQS3, 2), 2);
dEvS3_GeAz = cross(repmat(dEwL_GeAz,[1,1,nb]),rQS3,2) + cross(repmat(EwL_GeAz,[1,1,nb]), HvS3 + cross(repmat(EwL,[1,1,nb]), rQS3, 2), 2);
dEvS3_DrTr = cross(repmat(dEwL_DrTr,[1,1,nb]),rQS3,2) + cross(repmat(EwL_DrTr,[1,1,nb]), HvS3 + cross(repmat(EwL,[1,1,nb]), rQS3, 2), 2);
dEvS3_B3F1 = -coprod((s11_B3*qd_B3F1 + s12_B3*qd_B3F2 + s13_B3*qd_B3E1),J3(3,:)) + cross(repmat(EwL,1,1,nb),EvS3_B3F1,2);
dEvS3_B3E1 = -coprod((s33_B3*qd_B3E1 + s23_B3*qd_B3F2 + s13_B3*qd_B3F1),J3(3,:)) + cross(repmat(EwL,1,1,nb),EvS3_B3E1,2);
dEvS3_B3F2 = -coprod((s22_B3*qd_B3F2 + s12_B3*qd_B3F1 + s23_B3*qd_B3E1),J3(3,:)) + cross(repmat(EwL,1,1,nb),EvS3_B3F2,2);

dEvS3 = dEvS3_R*qd_R + dEvS3_P*qd_P + dEvS3_Y*qd_Y + dEvS3_TFA1*qd_TFA1 + dEvS3_TSS1*qd_TSS1 + dEvS3_TFA2*qd_TFA2 + dEvS3_TSS2*qd_TSS2 ...
      + dEvS3_yaw*qd_yaw + dEvS3_GeAz*qd_GeAz + dEvS3_DrTr*qd_DrTr + dEvS3_B3F1*qd_B3F1 + dEvS3_B3E1*qd_B3E1 + dEvS3_B3F2*qd_B3F2;

% end of Linear accelerations

%% Inertial Matrices
I_X = PtfmRIner*(A(1,:)'*A(1,:)) + PtfmYIner*(A(2,:)'*A(2,:)) + PtfmPIner*(A(3,:)'*A(3,:));
I_N = (NacYIner-mN*(NacCMxn^2+NacCMyn^2))*(D(2,:)'*D(2,:));
I_H = HubIner*(E(1,:)'*E(1,:));
I_G = GenIner*(C(1,:)'*C(1,:));

%% Partial forces and moments at Hub

% ==== Blade 1 Partial Forces ====
FB1_Sg   = - trapz(BldSec, coprod(mB1,EvS1_Sg), 3);
FB1_Sw   = - trapz(BldSec, coprod(mB1,EvS1_Sw), 3);
FB1_Hv   = - trapz(BldSec, coprod(mB1,EvS1_Hv), 3);
FB1_R    = - trapz(BldSec, coprod(mB1,EvS1_R), 3);
FB1_P    = - trapz(BldSec, coprod(mB1,EvS1_P), 3);
FB1_Y    = - trapz(BldSec, coprod(mB1,EvS1_Y), 3);
FB1_TFA1 = - trapz(BldSec, coprod(mB1,EvS1_TFA1), 3);
FB1_TSS1 = - trapz(BldSec, coprod(mB1,EvS1_TSS1), 3);
FB1_TFA2 = - trapz(BldSec, coprod(mB1,EvS1_TFA2), 3);
FB1_TSS2 = - trapz(BldSec, coprod(mB1,EvS1_TSS2), 3);
FB1_yaw  = - trapz(BldSec, coprod(mB1,EvS1_yaw), 3);
FB1_GeAz = - trapz(BldSec, coprod(mB1,EvS1_GeAz), 3);
FB1_DrTr = - trapz(BldSec, coprod(mB1,EvS1_DrTr), 3);
FB1_B1F1 = - trapz(BldSec, coprod(mB1,EvS1_B1F1), 3);
FB1_B1E1 = - trapz(BldSec, coprod(mB1,EvS1_B1E1), 3);
FB1_B1F2 = - trapz(BldSec, coprod(mB1,EvS1_B1F2), 3);

FB1_t    = trapz(BldSec, F_AeroB1 - coprod(mB1,(g*repmat(Z(2,:),[1,1,nb]) + dEvS1)), 3);

% ==== Blade 3 Partial Forces ====
FB2_Sg   = - trapz(BldSec, coprod(mB2,EvS2_Sg), 3);
FB2_Sw   = - trapz(BldSec, coprod(mB2,EvS2_Sw), 3);
FB2_Hv   = - trapz(BldSec, coprod(mB2,EvS2_Hv), 3);
FB2_R    = - trapz(BldSec, coprod(mB2,EvS2_R), 3);
FB2_P    = - trapz(BldSec, coprod(mB2,EvS2_P), 3);
FB2_Y    = - trapz(BldSec, coprod(mB2,EvS2_Y), 3);
FB2_TFA1 = - trapz(BldSec, coprod(mB2,EvS2_TFA1), 3);
FB2_TSS1 = - trapz(BldSec, coprod(mB2,EvS2_TSS1), 3);
FB2_TFA2 = - trapz(BldSec, coprod(mB2,EvS2_TFA2), 3);
FB2_TSS2 = - trapz(BldSec, coprod(mB2,EvS2_TSS2), 3);
FB2_yaw  = - trapz(BldSec, coprod(mB2,EvS2_yaw), 3);
FB2_GeAz = - trapz(BldSec, coprod(mB2,EvS2_GeAz), 3);
FB2_DrTr = - trapz(BldSec, coprod(mB2,EvS2_DrTr), 3);
FB2_B2F1 = - trapz(BldSec, coprod(mB2,EvS2_B2F1), 3);
FB2_B2E1 = - trapz(BldSec, coprod(mB2,EvS2_B2E1), 3);
FB2_B2F2 = - trapz(BldSec, coprod(mB2,EvS2_B2F2), 3);

FB2_t    = trapz(BldSec, F_AeroB2 - coprod(mB2,(g*repmat(Z(2,:),[1,1,nb]) + dEvS2)), 3);

% ==== Blade 3 Partial Forces ====
FB3_Sg   = - trapz(BldSec, coprod(mB2,EvS3_Sg), 3);
FB3_Sw   = - trapz(BldSec, coprod(mB2,EvS3_Sw), 3);
FB3_Hv   = - trapz(BldSec, coprod(mB2,EvS3_Hv), 3);
FB3_R    = - trapz(BldSec, coprod(mB2,EvS3_R), 3);
FB3_P    = - trapz(BldSec, coprod(mB2,EvS3_P), 3);
FB3_Y    = - trapz(BldSec, coprod(mB2,EvS3_Y), 3);
FB3_TFA1 = - trapz(BldSec, coprod(mB2,EvS3_TFA1), 3);
FB3_TSS1 = - trapz(BldSec, coprod(mB2,EvS3_TSS1), 3);
FB3_TFA2 = - trapz(BldSec, coprod(mB2,EvS3_TFA2), 3);
FB3_TSS2 = - trapz(BldSec, coprod(mB2,EvS3_TSS2), 3);
FB3_yaw  = - trapz(BldSec, coprod(mB2,EvS3_yaw), 3);
FB3_GeAz = - trapz(BldSec, coprod(mB2,EvS3_GeAz), 3);
FB3_DrTr = - trapz(BldSec, coprod(mB2,EvS3_DrTr), 3);
FB3_B3F1 = - trapz(BldSec, coprod(mB2,EvS3_B3F1), 3);
FB3_B3E1 = - trapz(BldSec, coprod(mB2,EvS3_B3E1), 3);
FB3_B3F2 = - trapz(BldSec, coprod(mB2,EvS3_B3F2), 3);

FB3_t    = trapz(BldSec, F_AeroB3 - coprod(mB2,(g*repmat(Z(2,:),[1,1,nb]) + dEvS3)), 3);

% ==== Blade 1 Partial Moments ====
% MB1_Sg   = - trapz(BldSec, cross(rQS1 - rQS10, mB1.*EvS1_Sg,   2), 3);
% MB1_Sw   = - trapz(BldSec, cross(rQS1 - rQS10, mB1.*EvS1_Sw,   2), 3);
% MB1_Hv   = - trapz(BldSec, cross(rQS1 - rQS10, mB1.*EvS1_Hv,   2), 3);
MB1_R    = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_R) ,   2), 3);
MB1_P    = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_P) ,   2), 3);
MB1_Y    = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_Y) ,   2), 3);
MB1_TFA1 = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_TFA1), 2), 3);
MB1_TSS1 = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_TSS1), 2), 3);
MB1_TFA2 = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_TFA2), 2), 3);
MB1_TSS2 = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_TSS2), 2), 3);
MB1_yaw  = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_yaw) , 2), 3);
MB1_GeAz = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_GeAz), 2), 3);
MB1_DrTr = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_DrTr), 2), 3);
MB1_B1F1 = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_B1F1), 2), 3);
MB1_B1E1 = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_B1E1), 2), 3);
MB1_B1F2 = - trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), coprod(mB1,EvS1_B1F2), 2), 3);

MB1_t    = trapz(BldSec, cross(rQS1 - repmat(rQS10,[1,1,nb]), F_AeroB1 - coprod(mB1,(g*repmat(Z(2,:),[1,1,nb]) + dEvS1)), 2), 3) + trapz(BldSec, M_AeroB1, 3);

% ==== Blade 2 Partial Moments ====
% MB2_Sg   = - trapz(BldSec, cross(rQS2 - rQS20, mB2.*EvS2_Sg,   2), 3);
% MB2_Sw   = - trapz(BldSec, cross(rQS2 - rQS20, mB2.*EvS2_Sw,   2), 3);
% MB2_Hv   = - trapz(BldSec, cross(rQS2 - rQS20, mB2.*EvS2_Hv,   2), 3);
MB2_R    = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_R),    2), 3);
MB2_P    = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_P),    2), 3);
MB2_Y    = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_Y),    2), 3);
MB2_TFA1 = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_TFA1), 2), 3);
MB2_TSS1 = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_TSS1), 2), 3);
MB2_TFA2 = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_TFA2), 2), 3);
MB2_TSS2 = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_TSS2), 2), 3);
MB2_yaw  = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_yaw) , 2), 3);
MB2_GeAz = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_GeAz), 2), 3);
MB2_DrTr = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_DrTr), 2), 3);
MB2_B2F1 = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_B2F1), 2), 3);
MB2_B2E1 = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_B2E1), 2), 3);
MB2_B2F2 = - trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), coprod(mB2,EvS2_B2F2), 2), 3);

MB2_t    = trapz(BldSec, cross(rQS2 - repmat(rQS20,[1,1,nb]), F_AeroB2 - coprod(mB2,(g*repmat(Z(2,:),[1,1,nb]) + dEvS2)), 2), 3) + trapz(BldSec, M_AeroB2, 3);

% ==== Blade 3 Partial Moments ====
% MB3_Sg   = - trapz(BldSec, cross(rQS3 - rQS30, mB3.*EvS3_Sg,   2), 3);
% MB3_Sw   = - trapz(BldSec, cross(rQS3 - rQS30, mB3.*EvS3_Sw,   2), 3);
% MB3_Hv   = - trapz(BldSec, cross(rQS3 - rQS30, mB3.*EvS3_Hv,   2), 3);
MB3_R    = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_R),    2), 3);
MB3_P    = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_P),    2), 3);
MB3_Y    = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_Y),    2), 3);
MB3_TFA1 = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_TFA1), 2), 3);
MB3_TSS1 = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_TSS1), 2), 3);
MB3_TFA2 = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_TFA2), 2), 3);
MB3_TSS2 = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_TSS2), 2), 3);
MB3_yaw  = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_yaw),  2), 3);
MB3_GeAz = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_GeAz), 2), 3);
MB3_DrTr = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_DrTr), 2), 3);
MB3_B3F1 = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_B3F1), 2), 3);
MB3_B3E1 = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_B3E1), 2), 3);
MB3_B3F2 = - trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), coprod(mB3,EvS3_B3F2), 2), 3);

MB3_t    = trapz(BldSec, cross(rQS3 - repmat(rQS30,[1,1,nb]), F_AeroB3 - coprod(mB3,(g*repmat(Z(2,:),[1,1,nb]) + dEvS3)), 2), 3) + trapz(BldSec, M_AeroB3, 3);

%% Partial forces and moments at Rotor
% Partial forces
FRotor_Sg   = FB1_Sg   + FB2_Sg   + FB3_Sg   - mH*EvC_Sg;
FRotor_Sw   = FB1_Sw   + FB2_Sw   + FB3_Sw   - mH*EvC_Sw;
FRotor_Hv   = FB1_Hv   + FB2_Hv   + FB3_Hv   - mH*EvC_Hv;
FRotor_R    = FB1_R    + FB2_R    + FB3_R    - mH*EvC_R;
FRotor_P    = FB1_P    + FB2_P    + FB3_P    - mH*EvC_P;
FRotor_Y    = FB1_Y    + FB2_Y    + FB3_Y    - mH*EvC_Y;
FRotor_TFA1 = FB1_TFA1 + FB2_TFA1 + FB3_TFA1 - mH*EvC_TFA1;
FRotor_TSS1 = FB1_TSS1 + FB2_TSS1 + FB3_TSS1 - mH*EvC_TSS1;
FRotor_TFA2 = FB1_TFA2 + FB2_TFA2 + FB3_TFA2 - mH*EvC_TFA2;
FRotor_TSS2 = FB1_TSS2 + FB2_TSS2 + FB3_TSS2 - mH*EvC_TSS2;
FRotor_yaw  = FB1_yaw  + FB2_yaw  + FB3_yaw  - mH*EvC_yaw;
FRotor_GeAz = FB1_GeAz + FB2_GeAz + FB3_GeAz - mH*EvC_GeAz;
FRotor_DrTr = FB1_DrTr + FB2_DrTr + FB3_DrTr - mH*EvC_DrTr;
FRotor_B1F1 = FB1_B1F1;
FRotor_B1E1 = FB1_B1E1;
FRotor_B1F2 = FB1_B1F2;
FRotor_B2F1 = FB2_B2F1;
FRotor_B2E1 = FB2_B2E1;
FRotor_B2F2 = FB2_B2F2;
FRotor_B3F1 = FB3_B3F1;
FRotor_B3E1 = FB3_B3E1;
FRotor_B3F2 = FB3_B3F2;

FRotor_t    = FB1_t + FB2_t + FB3_t - mH*(g*Z(2,:) + dEvC);

% Partial moments
% MRotor_Sg   = MB1_Sg   + MB2_Sg   + MB3_Sg    + cross(rQS10, FB1_Sg)   + cross(rQS20, FB2_Sg)   + cross(rQS30, FB3_Sg)   - mH*cross(rQC,EvC_Sg);
% MRotor_Sw   = MB1_Sw   + MB2_Sw   + MB3_Sw    + cross(rQS10, FB1_Sw)   + cross(rQS20, FB2_Sw)   + cross(rQS30, FB3_Sw)   - mH*cross(rQC,EvC_Sw);
% MRotor_Hv   = MB1_Hv   + MB2_Hv   + MB3_Hv    + cross(rQS10, FB1_Hv)   + cross(rQS20, FB2_Hv)   + cross(rQS30, FB3_Hv)   - mH*cross(rQC,EvC_Hv);
MRotor_R    = MB1_R    + MB2_R    + MB3_R     + cross(rQS10, FB1_R )   + cross(rQS20, FB2_R )   + cross(rQS30, FB3_R )   - mH*cross(rQC,EvC_R )   - (I_H*EwL_R' )';
MRotor_P    = MB1_P    + MB2_P    + MB3_P     + cross(rQS10, FB1_P )   + cross(rQS20, FB2_P )   + cross(rQS30, FB3_P )   - mH*cross(rQC,EvC_P )   - (I_H*EwL_P' )';
MRotor_Y    = MB1_Y    + MB2_Y    + MB3_Y     + cross(rQS10, FB1_Y )   + cross(rQS20, FB2_Y )   + cross(rQS30, FB3_Y )   - mH*cross(rQC,EvC_Y )   - (I_H*EwL_Y' )';
MRotor_TFA1 = MB1_TFA1 + MB2_TFA1 + MB3_TFA1  + cross(rQS10, FB1_TFA1) + cross(rQS20, FB2_TFA1) + cross(rQS30, FB3_TFA1) - mH*cross(rQC,EvC_TFA1) - (I_H*EwL_TFA1')';
MRotor_TSS1 = MB1_TSS1 + MB2_TSS1 + MB3_TSS1  + cross(rQS10, FB1_TSS1) + cross(rQS20, FB2_TSS1) + cross(rQS30, FB3_TSS1) - mH*cross(rQC,EvC_TSS1) - (I_H*EwL_TSS1')';
MRotor_TFA2 = MB1_TFA2 + MB2_TFA2 + MB3_TFA2  + cross(rQS10, FB1_TFA2) + cross(rQS20, FB2_TFA2) + cross(rQS30, FB3_TFA2) - mH*cross(rQC,EvC_TFA2) - (I_H*EwL_TFA2')';
MRotor_TSS2 = MB1_TSS2 + MB2_TSS2 + MB3_TSS2  + cross(rQS10, FB1_TSS2) + cross(rQS20, FB2_TSS2) + cross(rQS30, FB3_TSS2) - mH*cross(rQC,EvC_TSS2) - (I_H*EwL_TSS2')';
MRotor_yaw  = MB1_yaw  + MB2_yaw  + MB3_yaw   + cross(rQS10, FB1_yaw)  + cross(rQS20, FB2_yaw)  + cross(rQS30, FB3_yaw)  - mH*cross(rQC,EvC_yaw)  - (I_H*EwL_yaw')';
MRotor_GeAz = MB1_GeAz + MB2_GeAz + MB3_GeAz  + cross(rQS10, FB1_GeAz) + cross(rQS20, FB2_GeAz) + cross(rQS30, FB3_GeAz) - mH*cross(rQC,EvC_GeAz) - (I_H*EwL_GeAz')';
MRotor_DrTr = MB1_DrTr + MB2_DrTr + MB3_DrTr  + cross(rQS10, FB1_DrTr) + cross(rQS20, FB2_DrTr) + cross(rQS30, FB3_DrTr) - mH*cross(rQC,EvC_DrTr) - (I_H*EwL_DrTr')';
MRotor_B1F1 = MB1_B1F1 + cross(rQS10, FB1_B1F1);
MRotor_B1E1 = MB1_B1E1 + cross(rQS10, FB1_B1E1);
MRotor_B1F2 = MB1_B1F2 + cross(rQS10, FB1_B1F2);
MRotor_B2F1 = MB2_B2F1 + cross(rQS20, FB2_B2F1);
MRotor_B2E1 = MB2_B2E1 + cross(rQS20, FB2_B2E1);
MRotor_B2F2 = MB2_B2F2 + cross(rQS20, FB2_B2F2);
MRotor_B3F1 = MB3_B3F1 + cross(rQS30, FB3_B3F1);
MRotor_B3E1 = MB3_B3E1 + cross(rQS30, FB3_B3E1);
MRotor_B3F2 = MB3_B3F2 + cross(rQS30, FB3_B3F2);

MRotor_t    = MB1_t + MB2_t + MB3_t  + cross(rQS10, FB1_t) + cross(rQS20, FB2_t) + cross(rQS30, FB3_t) - mH*cross(rQC, dEvC + g*Z(2,:)) - (I_H*dEwL')' - cross(EwL, (I_H*EwL')');

%% Partial forces and moments at tower top
% Patial forces
FO_Sg   = FRotor_Sg   - mN*EvU_Sg;
FO_Sw   = FRotor_Sw   - mN*EvU_Sw;
FO_Hv   = FRotor_Hv   - mN*EvU_Hv;
FO_R    = FRotor_R    - mN*EvU_R ;
FO_P    = FRotor_P    - mN*EvU_P;
FO_Y    = FRotor_Y    - mN*EvU_Y;
FO_TFA1 = FRotor_TFA1 - mN*EvU_TFA1;
FO_TSS1 = FRotor_TSS1 - mN*EvU_TSS1;
FO_TFA2 = FRotor_TFA2 - mN*EvU_TFA2;
FO_TSS2 = FRotor_TSS2 - mN*EvU_TSS2;
FO_yaw  = FRotor_yaw  - mN*EvU_yaw;
FO_GeAz = FRotor_GeAz;
FO_DrTr = FRotor_DrTr;
FO_B1F1 = FRotor_B1F1;
FO_B1E1 = FRotor_B1E1;
FO_B1F2 = FRotor_B1F2;
FO_B2F1 = FRotor_B2F1;
FO_B2E1 = FRotor_B2E1;
FO_B2F2 = FRotor_B2F2;
FO_B3F1 = FRotor_B3F1;
FO_B3E1 = FRotor_B3E1;
FO_B3F2 = FRotor_B3F2;

FO_t    = FRotor_t - mN*(dEvU + g*Z(2,:));

% Partial moments
% MO_Sg   = MRotor_Sg   + cross(rOQ, FRotor_Sg)	                   - mN*cross(rOU,EvU_Sg);
% MO_Sw   = MRotor_Sw   + cross(rOQ, FRotor_Sw)                      - mN*cross(rOU,EvU_Sw);
% MO_Hv   = MRotor_Hv   + cross(rOQ, FRotor_Hv)                      - mN*cross(rOU,EvU_Hv);
MO_R    = MRotor_R    + cross(rOQ, FRotor_R )   - (I_G*EwG_R' )'   - mN*cross(rOU,EvU_R )   - (I_N*EwN_R' )';
MO_P    = MRotor_P    + cross(rOQ, FRotor_P )   - (I_G*EwG_P' )'   - mN*cross(rOU,EvU_P )   - (I_N*EwN_P' )';
MO_Y    = MRotor_Y    + cross(rOQ, FRotor_Y )   - (I_G*EwG_Y' )'   - mN*cross(rOU,EvU_Y )   - (I_N*EwN_Y' )';
MO_TFA1 = MRotor_TFA1 + cross(rOQ, FRotor_TFA1) - (I_G*EwG_TFA1')' - mN*cross(rOU,EvU_TFA1) - (I_N*EwN_TFA1')';
MO_TSS1 = MRotor_TSS1 + cross(rOQ, FRotor_TSS1) - (I_G*EwG_TSS1')' - mN*cross(rOU,EvU_TSS1) - (I_N*EwN_TSS1')';
MO_TFA2 = MRotor_TFA2 + cross(rOQ, FRotor_TFA2) - (I_G*EwG_TFA2')' - mN*cross(rOU,EvU_TFA2) - (I_N*EwN_TFA2')';
MO_TSS2 = MRotor_TSS2 + cross(rOQ, FRotor_TSS2) - (I_G*EwG_TSS2')' - mN*cross(rOU,EvU_TSS2) - (I_N*EwN_TSS2')';
MO_yaw  = MRotor_yaw  + cross(rOQ, FRotor_yaw)  - (I_G*EwG_yaw')'  - mN*cross(rOU,EvU_yaw)  - (I_N*EwN_yaw')';
MO_GeAz = MRotor_GeAz + cross(rOQ, FRotor_GeAz) - (I_G*EwG_GeAz')';
MO_DrTr = MRotor_DrTr + cross(rOQ, FRotor_DrTr);
MO_B1F1 = MRotor_B1F1 + cross(rOQ, FRotor_B1F1);
MO_B1E1 = MRotor_B1E1 + cross(rOQ, FRotor_B1E1);
MO_B1F2 = MRotor_B1F2 + cross(rOQ, FRotor_B1F2);
MO_B2F1 = MRotor_B2F1 + cross(rOQ, FRotor_B2F1);
MO_B2E1 = MRotor_B2E1 + cross(rOQ, FRotor_B2E1);
MO_B2F2 = MRotor_B2F2 + cross(rOQ, FRotor_B2F2);
MO_B3F1 = MRotor_B3F1 + cross(rOQ, FRotor_B3F1);
MO_B3E1 = MRotor_B3E1 + cross(rOQ, FRotor_B3E1);
MO_B3F2 = MRotor_B3F2 + cross(rOQ, FRotor_B3F2);

MO_t    = MRotor_t + cross(rOQ, FRotor_t) - (I_G*dEwG')' - cross(EwG, (I_G*EwG')') - mN*cross(rOU, dEvU + g*Z(2,:)) - (I_N*dEwN')' - cross(EwN, (I_N*EwN')');

%% Partial forces and moments at tower base
% Partial forces
FT_Sg   = FO_Sg   - YawBrMass*EvO_Sg   - trapz(TwrSec, coprod(mT,EvF_Sg), 3);
FT_Sw   = FO_Sw   - YawBrMass*EvO_Sw   - trapz(TwrSec, coprod(mT,EvF_Sw), 3);
FT_Hv   = FO_Hv   - YawBrMass*EvO_Hv   - trapz(TwrSec, coprod(mT,EvF_Hv), 3);
FT_R    = FO_R    - YawBrMass*EvO_R    - trapz(TwrSec, coprod(mT,EvF_R ), 3);
FT_P    = FO_P    - YawBrMass*EvO_P    - trapz(TwrSec, coprod(mT,EvF_P ), 3);
FT_Y    = FO_Y    - YawBrMass*EvO_Y    - trapz(TwrSec, coprod(mT,EvF_Y ), 3);
FT_TFA1 = FO_TFA1 - YawBrMass*EvO_TFA1 - trapz(TwrSec, coprod(mT,EvF_TFA1), 3);
FT_TSS1 = FO_TSS1 - YawBrMass*EvO_TSS1 - trapz(TwrSec, coprod(mT,EvF_TSS1), 3);
FT_TFA2 = FO_TFA2 - YawBrMass*EvO_TFA2 - trapz(TwrSec, coprod(mT,EvF_TFA2), 3);
FT_TSS2 = FO_TSS2 - YawBrMass*EvO_TSS2 - trapz(TwrSec, coprod(mT,EvF_TSS2), 3);
FT_yaw  = FO_yaw;
FT_GeAz = FO_GeAz;
FT_DrTr = FO_DrTr;
FT_B1F1 = FO_B1F1;
FT_B1E1 = FO_B1E1;
FT_B1F2 = FO_B1F2;
FT_B2F1 = FO_B2F1;
FT_B2E1 = FO_B2E1;
FT_B2F2 = FO_B2F2;
FT_B3F1 = FO_B3F1;
FT_B3E1 = FO_B3E1;
FT_B3F2 = FO_B3F2;

FT_t    = FO_t  - YawBrMass*(dEvO + g*Z(2,:)) + trapz(TwrSec, F_AeroT - coprod(mT,(dEvF + g*repmat(Z(2,:),[1,1,nt]))), 3);

% Partial moments
% MT_Sg   = MO_Sg   + cross(rZO - rZT0, FO_Sg   - YawBrMass*EvO_Sg)   - trapz(TwrSec, cross(rZT - rZT0, mT.*EvF_Sg, 2), 3);
% MT_Sw   = MO_Sw   + cross(rZO - rZT0, FO_Sw   - YawBrMass*EvO_Sw)   - trapz(TwrSec, cross(rZT - rZT0, mT.*EvF_Sw, 2), 3);
% MT_Hv   = MO_Hv   + cross(rZO - rZT0, FO_Hv   - YawBrMass*EvO_Hv)   - trapz(TwrSec, cross(rZT - rZT0, mT.*EvF_Hv, 2), 3);
MT_R    = MO_R    + cross(rZO - rZT0, FO_R    - YawBrMass*EvO_R )   - trapz(TwrSec, cross(rZT - repmat(rZT0,[1,1,nt]), coprod(mT,EvF_R) , 2), 3);
MT_P    = MO_P    + cross(rZO - rZT0, FO_P    - YawBrMass*EvO_P )   - trapz(TwrSec, cross(rZT - repmat(rZT0,[1,1,nt]), coprod(mT,EvF_P) , 2), 3);
MT_Y    = MO_Y    + cross(rZO - rZT0, FO_Y    - YawBrMass*EvO_Y )   - trapz(TwrSec, cross(rZT - repmat(rZT0,[1,1,nt]), coprod(mT,EvF_Y) , 2), 3);
MT_TFA1 = MO_TFA1 + cross(rZO - rZT0, FO_TFA1 - YawBrMass*EvO_TFA1) - trapz(TwrSec, cross(rZT - repmat(rZT0,[1,1,nt]), coprod(mT,EvF_TFA1), 2), 3);
MT_TSS1 = MO_TSS1 + cross(rZO - rZT0, FO_TSS1 - YawBrMass*EvO_TSS1) - trapz(TwrSec, cross(rZT - repmat(rZT0,[1,1,nt]), coprod(mT,EvF_TSS1), 2), 3);
MT_TFA2 = MO_TFA2 + cross(rZO - rZT0, FO_TFA2 - YawBrMass*EvO_TFA2) - trapz(TwrSec, cross(rZT - repmat(rZT0,[1,1,nt]), coprod(mT,EvF_TFA2), 2), 3);
MT_TSS2 = MO_TSS2 + cross(rZO - rZT0, FO_TSS2 - YawBrMass*EvO_TSS2) - trapz(TwrSec, cross(rZT - repmat(rZT0,[1,1,nt]), coprod(mT,EvF_TSS2), 2), 3);
MT_yaw  = MO_yaw  + cross(rZO - rZT0, FO_yaw);
MT_GeAz = MO_GeAz + cross(rZO - rZT0, FO_GeAz);
MT_DrTr = MO_DrTr + cross(rZO - rZT0, FO_DrTr);
MT_B1F1 = MO_B1F1 + cross(rZO - rZT0, FO_B1F1);
MT_B1E1 = MO_B1E1 + cross(rZO - rZT0, FO_B1E1);
MT_B1F2 = MO_B1F2 + cross(rZO - rZT0, FO_B1F2);
MT_B2F1 = MO_B2F1 + cross(rZO - rZT0, FO_B2F1);
MT_B2E1 = MO_B2E1 + cross(rZO - rZT0, FO_B2E1);
MT_B2F2 = MO_B2F2 + cross(rZO - rZT0, FO_B2F2);
MT_B3F1 = MO_B3F1 + cross(rZO - rZT0, FO_B3F1);
MT_B3E1 = MO_B3E1 + cross(rZO - rZT0, FO_B3E1);
MT_B3F2 = MO_B3F2 + cross(rZO - rZT0, FO_B3F2);

MT_t    = MO_t + cross(rZO - rZT0, FO_t - YawBrMass*(dEvO + g*Z(2,:))) + trapz(TwrSec, cross(rZT - repmat(rZT0,[1,1,nt]), F_AeroT - coprod(mT,(dEvF + g*repmat(Z(2,:),[1,1,nt]))), 2), 3);

%% Partial forces and moments at platform CM
% Partial forces
FAll_Sg   = FT_Sg - mX*EvY_Sg;
FAll_Sw   = FT_Sw - mX*EvY_Sw;
FAll_Hv   = FT_Hv - mX*EvY_Hv;
FAll_R    = FT_R  - mX*EvY_R ;
FAll_P    = FT_P  - mX*EvY_P ;
FAll_Y    = FT_Y  - mX*EvY_Y ;
FAll_TFA1 = FT_TFA1;
FAll_TSS1 = FT_TSS1;
FAll_TFA2 = FT_TFA2;
FAll_TSS2 = FT_TSS2;
FAll_yaw  = FT_yaw;
FAll_GeAz = FT_GeAz;
FAll_DrTr = FT_DrTr;
FAll_B1F1 = FT_B1F1;
FAll_B1E1 = FT_B1E1;
FAll_B1F2 = FT_B1F2;
FAll_B2F1 = FT_B2F1;
FAll_B2E1 = FT_B2E1;
FAll_B2F2 = FT_B2F2;
FAll_B3F1 = FT_B3F1;
FAll_B3E1 = FT_B3E1;
FAll_B3F2 = FT_B3F2;

FAll_t    = FT_t + F_Hydro - mX*(dEvY + g*Z(2,:));

% Partial moments
% MAll_Sg   = MT_Sg   + cross(rZT0, FT_Sg) - mX*cross(rZY, EvY_Sg);
% MAll_Sw   = MT_Sw   + cross(rZT0, FT_Sw) - mX*cross(rZY, EvY_Sw);
% MAll_Hv   = MT_Hv   + cross(rZT0, FT_Hv) - mX*cross(rZY, EvY_Hv);
MAll_R    = MT_R    + cross(rZT0, FT_R ) - mX*cross(rZY, EvY_R ) - (I_X*EwX_R' )';
MAll_P    = MT_P    + cross(rZT0, FT_P ) - mX*cross(rZY, EvY_P ) - (I_X*EwX_P' )';
MAll_Y    = MT_Y    + cross(rZT0, FT_Y ) - mX*cross(rZY, EvY_Y ) - (I_X*EwX_Y' )';
MAll_TFA1 = MT_TFA1 + cross(rZT0, FT_TFA1); 
MAll_TSS1 = MT_TSS1 + cross(rZT0, FT_TSS1);
MAll_TFA2 = MT_TFA2 + cross(rZT0, FT_TFA2);
MAll_TSS2 = MT_TSS2 + cross(rZT0, FT_TSS2);
MAll_yaw  = MT_yaw  + cross(rZT0, FT_yaw);
MAll_GeAz = MT_GeAz + cross(rZT0, FT_GeAz);
MAll_DrTr = MT_DrTr + cross(rZT0, FT_DrTr);
MAll_B1F1 = MT_B1F1 + cross(rZT0, FT_B1F1);
MAll_B1E1 = MT_B1E1 + cross(rZT0, FT_B1E1);
MAll_B1F2 = MT_B1F2 + cross(rZT0, FT_B1F2);
MAll_B2F1 = MT_B2F1 + cross(rZT0, FT_B2F1);
MAll_B2E1 = MT_B2E1 + cross(rZT0, FT_B2E1);
MAll_B2F2 = MT_B2F2 + cross(rZT0, FT_B2F2);
MAll_B3F1 = MT_B3F1 + cross(rZT0, FT_B3F1);
MAll_B3E1 = MT_B3E1 + cross(rZT0, FT_B3E1);
MAll_B3F2 = MT_B3F2 + cross(rZT0, FT_B3F2);

MAll_t    = MT_t + cross(rZT0, FT_t) + M_Hydro - mX*cross(rZY, dEvY + g*Z(2,:)) - cross(EwX, (I_X*EwX')');

%% Initialize matrices
IM_nom   = zeros(22,22);
C_B1 = IM_nom; C_B2 = IM_nom; C_B3 = IM_nom; C_T = IM_nom;
C_com = IM_nom;

f_com = zeros(22,1);
f_B1 = f_com; f_ElasticB1 = f_com; f_DampingB1 = f_com;  f_B2 = f_com; f_ElasticB2 = f_com; f_DampingB2 = f_com; f_B3 = f_com; f_ElasticB3 = f_com; f_DampingB3 = f_com; 
f_Yaw = f_com; f_T = f_com; f_ElasticT = f_com; f_DampingT = f_com; f_G = f_com; f_ElasticDrive = f_com; f_DampDrive = f_com; f_Gen = f_com; f_Brake = f_com;

%% Common C terms

% r = 1, 2, 3
C_com(1,1)  = -dot(EvZ_Sg, FAll_Sg);
C_com(1,2)  = -dot(EvZ_Sg, FAll_Sw);
C_com(1,3)  = -dot(EvZ_Sg, FAll_Hv);
C_com(1,4)  = -dot(EvZ_Sg, FAll_R );
C_com(1,5)  = -dot(EvZ_Sg, FAll_P );
C_com(1,6)  = -dot(EvZ_Sg, FAll_Y );
C_com(1,7)  = -dot(EvZ_Sg, FAll_TFA1);
C_com(1,8)  = -dot(EvZ_Sg, FAll_TSS1);
C_com(1,9)  = -dot(EvZ_Sg, FAll_TFA2);
C_com(1,10) = -dot(EvZ_Sg, FAll_TSS2);
C_com(1,11) = -dot(EvZ_Sg, FAll_yaw );
C_com(1,12) = -dot(EvZ_Sg, FAll_GeAz);
C_com(1,13) = -dot(EvZ_Sg, FAll_DrTr);
C_com(1,14) = -dot(EvZ_Sg, FAll_B1F1);
C_com(1,15) = -dot(EvZ_Sg, FAll_B1E1);
C_com(1,16) = -dot(EvZ_Sg, FAll_B1F2);
C_com(1,17) = -dot(EvZ_Sg, FAll_B2F1);
C_com(1,18) = -dot(EvZ_Sg, FAll_B2E1);
C_com(1,19) = -dot(EvZ_Sg, FAll_B2F2);
C_com(1,20) = -dot(EvZ_Sg, FAll_B3F1);
C_com(1,21) = -dot(EvZ_Sg, FAll_B3E1);
C_com(1,22) = -dot(EvZ_Sg, FAll_B3F2);

C_com(2,2)  = -dot(EvZ_Sw, FAll_Sw);
C_com(2,3)  = -dot(EvZ_Sw, FAll_Hv);
C_com(2,4)  = -dot(EvZ_Sw, FAll_R );
C_com(2,5)  = -dot(EvZ_Sw, FAll_P );
C_com(2,6)  = -dot(EvZ_Sw, FAll_Y );
C_com(2,7)  = -dot(EvZ_Sw, FAll_TFA1);
C_com(2,8)  = -dot(EvZ_Sw, FAll_TSS1);
C_com(2,9)  = -dot(EvZ_Sw, FAll_TFA2);
C_com(2,10) = -dot(EvZ_Sw, FAll_TSS2);
C_com(2,11) = -dot(EvZ_Sw, FAll_yaw );
C_com(2,12) = -dot(EvZ_Sw, FAll_GeAz);
C_com(2,13) = -dot(EvZ_Sw, FAll_DrTr);
C_com(2,14) = -dot(EvZ_Sw, FAll_B1F1);
C_com(2,15) = -dot(EvZ_Sw, FAll_B1E1);
C_com(2,16) = -dot(EvZ_Sw, FAll_B1F2);
C_com(2,17) = -dot(EvZ_Sw, FAll_B2F1);
C_com(2,18) = -dot(EvZ_Sw, FAll_B2E1);
C_com(2,19) = -dot(EvZ_Sw, FAll_B2F2);
C_com(2,20) = -dot(EvZ_Sw, FAll_B3F1);
C_com(2,21) = -dot(EvZ_Sw, FAll_B3E1);
C_com(2,22) = -dot(EvZ_Sw, FAll_B3F2);

C_com(3,3)  = -dot(EvZ_Hv, FAll_Hv);
C_com(3,4)  = -dot(EvZ_Hv, FAll_R );
C_com(3,5)  = -dot(EvZ_Hv, FAll_P );
C_com(3,6)  = -dot(EvZ_Hv, FAll_Y );
C_com(3,7)  = -dot(EvZ_Hv, FAll_TFA1);
C_com(3,8)  = -dot(EvZ_Hv, FAll_TSS1);
C_com(3,9)  = -dot(EvZ_Hv, FAll_TFA2);
C_com(3,10) = -dot(EvZ_Hv, FAll_TSS2);
C_com(3,11) = -dot(EvZ_Hv, FAll_yaw );
C_com(3,12) = -dot(EvZ_Hv, FAll_GeAz);
C_com(3,13) = -dot(EvZ_Hv, FAll_DrTr);
C_com(3,14) = -dot(EvZ_Hv, FAll_B1F1);
C_com(3,15) = -dot(EvZ_Hv, FAll_B1E1);
C_com(3,16) = -dot(EvZ_Hv, FAll_B1F2);
C_com(3,17) = -dot(EvZ_Hv, FAll_B2F1);
C_com(3,18) = -dot(EvZ_Hv, FAll_B2E1);
C_com(3,19) = -dot(EvZ_Hv, FAll_B2F2);
C_com(3,20) = -dot(EvZ_Hv, FAll_B3F1);
C_com(3,21) = -dot(EvZ_Hv, FAll_B3E1);
C_com(3,22) = -dot(EvZ_Hv, FAll_B3F2);

% r = 4, 5 ,6
C_com(4,4)  = -dot(EwX_R, MAll_R );
C_com(4,5)  = -dot(EwX_R, MAll_P );
C_com(4,6)  = -dot(EwX_R, MAll_Y );
C_com(4,7)  = -dot(EwX_R, MAll_TFA1);
C_com(4,8)  = -dot(EwX_R, MAll_TSS1);
C_com(4,9)  = -dot(EwX_R, MAll_TFA2);
C_com(4,10) = -dot(EwX_R, MAll_TSS2);
C_com(4,11) = -dot(EwX_R, MAll_yaw );
C_com(4,12) = -dot(EwX_R, MAll_GeAz);
C_com(4,13) = -dot(EwX_R, MAll_DrTr);
C_com(4,14) = -dot(EwX_R, MAll_B1F1);
C_com(4,15) = -dot(EwX_R, MAll_B1E1);
C_com(4,16) = -dot(EwX_R, MAll_B1F2);
C_com(4,17) = -dot(EwX_R, MAll_B2F1);
C_com(4,18) = -dot(EwX_R, MAll_B2E1);
C_com(4,19) = -dot(EwX_R, MAll_B2F2);
C_com(4,20) = -dot(EwX_R, MAll_B3F1);
C_com(4,21) = -dot(EwX_R, MAll_B3E1);
C_com(4,22) = -dot(EwX_R, MAll_B3F2);

C_com(5,5)  = -dot(EwX_P, MAll_P );
C_com(5,6)  = -dot(EwX_P, MAll_Y );
C_com(5,7)  = -dot(EwX_P, MAll_TFA1);
C_com(5,8)  = -dot(EwX_P, MAll_TSS1);
C_com(5,9)  = -dot(EwX_P, MAll_TFA2);
C_com(5,10) = -dot(EwX_P, MAll_TSS2);
C_com(5,11) = -dot(EwX_P, MAll_yaw );
C_com(5,12) = -dot(EwX_P, MAll_GeAz);
C_com(5,13) = -dot(EwX_P, MAll_DrTr);
C_com(5,14) = -dot(EwX_P, MAll_B1F1);
C_com(5,15) = -dot(EwX_P, MAll_B1E1);
C_com(5,16) = -dot(EwX_P, MAll_B1F2);
C_com(5,17) = -dot(EwX_P, MAll_B2F1);
C_com(5,18) = -dot(EwX_P, MAll_B2E1);
C_com(5,19) = -dot(EwX_P, MAll_B2F2);
C_com(5,20) = -dot(EwX_P, MAll_B3F1);
C_com(5,21) = -dot(EwX_P, MAll_B3E1);
C_com(5,22) = -dot(EwX_P, MAll_B3F2);

C_com(6,6)  = -dot(EwX_Y, MAll_Y );
C_com(6,7)  = -dot(EwX_Y, MAll_TFA1);
C_com(6,8)  = -dot(EwX_Y, MAll_TSS1);
C_com(6,9)  = -dot(EwX_Y, MAll_TFA2);
C_com(6,10) = -dot(EwX_Y, MAll_TSS2);
C_com(6,11) = -dot(EwX_Y, MAll_yaw );
C_com(6,12) = -dot(EwX_Y, MAll_GeAz);
C_com(6,13) = -dot(EwX_Y, MAll_DrTr);
C_com(6,14) = -dot(EwX_Y, MAll_B1F1);
C_com(6,15) = -dot(EwX_Y, MAll_B1E1);
C_com(6,16) = -dot(EwX_Y, MAll_B1F2);
C_com(6,17) = -dot(EwX_Y, MAll_B2F1);
C_com(6,18) = -dot(EwX_Y, MAll_B2E1);
C_com(6,19) = -dot(EwX_Y, MAll_B2F2);
C_com(6,20) = -dot(EwX_Y, MAll_B3F1);
C_com(6,21) = -dot(EwX_Y, MAll_B3E1);
C_com(6,22) = -dot(EwX_Y, MAll_B3F2);

% r = 7, 8, 9, 10
C_com(7,7)  = - dot(EvO_TFA1,FO_TFA1) - dot(EwB_TFA1, MO_TFA1);
C_com(7,8)  = - dot(EvO_TFA1,FO_TSS1) - dot(EwB_TFA1, MO_TSS1);
C_com(7,9)  = - dot(EvO_TFA1,FO_TFA2) - dot(EwB_TFA1, MO_TFA2);
C_com(7,10) = - dot(EvO_TFA1,FO_TSS2) - dot(EwB_TFA1, MO_TSS2);
C_com(7,11) = - dot(EvO_TFA1,FO_yaw ) - dot(EwB_TFA1, MO_yaw );
C_com(7,12) = - dot(EvO_TFA1,FO_GeAz) - dot(EwB_TFA1, MO_GeAz);
C_com(7,13) = - dot(EvO_TFA1,FO_DrTr) - dot(EwB_TFA1, MO_DrTr);
C_com(7,14) = - dot(EvO_TFA1,FO_B1F1) - dot(EwB_TFA1, MO_B1F1);
C_com(7,15) = - dot(EvO_TFA1,FO_B1E1) - dot(EwB_TFA1, MO_B1E1);
C_com(7,16) = - dot(EvO_TFA1,FO_B1F2) - dot(EwB_TFA1, MO_B1F2);
C_com(7,17) = - dot(EvO_TFA1,FO_B2F1) - dot(EwB_TFA1, MO_B2F1);
C_com(7,18) = - dot(EvO_TFA1,FO_B2E1) - dot(EwB_TFA1, MO_B2E1);
C_com(7,19) = - dot(EvO_TFA1,FO_B2F2) - dot(EwB_TFA1, MO_B2F2);
C_com(7,20) = - dot(EvO_TFA1,FO_B3F1) - dot(EwB_TFA1, MO_B3F1);
C_com(7,21) = - dot(EvO_TFA1,FO_B3E1) - dot(EwB_TFA1, MO_B3E1);
C_com(7,22) = - dot(EvO_TFA1,FO_B3F2) - dot(EwB_TFA1, MO_B3F2);

C_com(8,8)  = - dot(EvO_TSS1,FO_TSS1) - dot(EwB_TSS1, MO_TSS1);
C_com(8,9)  = - dot(EvO_TSS1,FO_TFA2) - dot(EwB_TSS1, MO_TFA2);
C_com(8,10) = - dot(EvO_TSS1,FO_TSS2) - dot(EwB_TSS1, MO_TSS2);
C_com(8,11) = - dot(EvO_TSS1,FO_yaw ) - dot(EwB_TSS1, MO_yaw );
C_com(8,12) = - dot(EvO_TSS1,FO_GeAz) - dot(EwB_TSS1, MO_GeAz);
C_com(8,13) = - dot(EvO_TSS1,FO_DrTr) - dot(EwB_TSS1, MO_DrTr);
C_com(8,14) = - dot(EvO_TSS1,FO_B1F1) - dot(EwB_TSS1, MO_B1F1);
C_com(8,15) = - dot(EvO_TSS1,FO_B1E1) - dot(EwB_TSS1, MO_B1E1);
C_com(8,16) = - dot(EvO_TSS1,FO_B1F2) - dot(EwB_TSS1, MO_B1F2);
C_com(8,17) = - dot(EvO_TSS1,FO_B2F1) - dot(EwB_TSS1, MO_B2F1);
C_com(8,18) = - dot(EvO_TSS1,FO_B2E1) - dot(EwB_TSS1, MO_B2E1);
C_com(8,19) = - dot(EvO_TSS1,FO_B2F2) - dot(EwB_TSS1, MO_B2F2);
C_com(8,20) = - dot(EvO_TSS1,FO_B3F1) - dot(EwB_TSS1, MO_B3F1);
C_com(8,21) = - dot(EvO_TSS1,FO_B3E1) - dot(EwB_TSS1, MO_B3E1);
C_com(8,22) = - dot(EvO_TSS1,FO_B3F2) - dot(EwB_TSS1, MO_B3F2);

C_com(9,9)  = - dot(EvO_TFA2,FO_TFA2) - dot(EwB_TFA2, MO_TFA2);
C_com(9,10) = - dot(EvO_TFA2,FO_TSS2) - dot(EwB_TFA2, MO_TSS2);
C_com(9,11) = - dot(EvO_TFA2,FO_yaw ) - dot(EwB_TFA2, MO_yaw );
C_com(9,12) = - dot(EvO_TFA2,FO_GeAz) - dot(EwB_TFA2, MO_GeAz);
C_com(9,13) = - dot(EvO_TFA2,FO_DrTr) - dot(EwB_TFA2, MO_DrTr);
C_com(9,14) = - dot(EvO_TFA2,FO_B1F1) - dot(EwB_TFA2, MO_B1F1);
C_com(9,15) = - dot(EvO_TFA2,FO_B1E1) - dot(EwB_TFA2, MO_B1E1);
C_com(9,16) = - dot(EvO_TFA2,FO_B1F2) - dot(EwB_TFA2, MO_B1F2);
C_com(9,17) = - dot(EvO_TFA2,FO_B2F1) - dot(EwB_TFA2, MO_B2F1);
C_com(9,18) = - dot(EvO_TFA2,FO_B2E1) - dot(EwB_TFA2, MO_B2E1);
C_com(9,19) = - dot(EvO_TFA2,FO_B2F2) - dot(EwB_TFA2, MO_B2F2);
C_com(9,20) = - dot(EvO_TFA2,FO_B3F1) - dot(EwB_TFA2, MO_B3F1);
C_com(9,21) = - dot(EvO_TFA2,FO_B3E1) - dot(EwB_TFA2, MO_B3E1);
C_com(9,22) = - dot(EvO_TFA2,FO_B3F2) - dot(EwB_TFA2, MO_B3F2);

C_com(10,10) = - dot(EvO_TSS2,FO_TSS2) - dot(EwB_TSS2, MO_TSS2);
C_com(10,11) = - dot(EvO_TSS2,FO_yaw ) - dot(EwB_TSS2, MO_yaw );
C_com(10,12) = - dot(EvO_TSS2,FO_GeAz) - dot(EwB_TSS2, MO_GeAz);
C_com(10,13) = - dot(EvO_TSS2,FO_DrTr) - dot(EwB_TSS2, MO_DrTr);
C_com(10,14) = - dot(EvO_TSS2,FO_B1F1) - dot(EwB_TSS2, MO_B1F1);
C_com(10,15) = - dot(EvO_TSS2,FO_B1E1) - dot(EwB_TSS2, MO_B1E1);
C_com(10,16) = - dot(EvO_TSS2,FO_B1F2) - dot(EwB_TSS2, MO_B1F2);
C_com(10,17) = - dot(EvO_TSS2,FO_B2F1) - dot(EwB_TSS2, MO_B2F1);
C_com(10,18) = - dot(EvO_TSS2,FO_B2E1) - dot(EwB_TSS2, MO_B2E1);
C_com(10,19) = - dot(EvO_TSS2,FO_B2F2) - dot(EwB_TSS2, MO_B2F2);
C_com(10,20) = - dot(EvO_TSS2,FO_B3F1) - dot(EwB_TSS2, MO_B3F1);
C_com(10,21) = - dot(EvO_TSS2,FO_B3E1) - dot(EwB_TSS2, MO_B3E1);
C_com(10,22) = - dot(EvO_TSS2,FO_B3F2) - dot(EwB_TSS2, MO_B3F2);

% r = 11
C_com(11,11) = -dot(EwN_yaw, MO_yaw); 
C_com(11,12) = -dot(EwN_yaw, MO_GeAz); 
C_com(11,13) = -dot(EwN_yaw, MO_DrTr); 
C_com(11,14) = -dot(EwN_yaw, MO_B1F1); 
C_com(11,15) = -dot(EwN_yaw, MO_B1E1);
C_com(11,16) = -dot(EwN_yaw, MO_B1F2); 
C_com(11,17) = -dot(EwN_yaw, MO_B2F1); 
C_com(11,18) = -dot(EwN_yaw, MO_B2E1); 
C_com(11,19) = -dot(EwN_yaw, MO_B2F2); 
C_com(11,20) = -dot(EwN_yaw, MO_B3F1);
C_com(11,21) = -dot(EwN_yaw, MO_B3E1); 
C_com(11,22) = -dot(EwN_yaw, MO_B3F2);

% r = 12   
C_com(12,12) = -dot(EwL_GeAz, MRotor_GeAz) + dot(EwG_GeAz,I_G*EwG_GeAz'); 
C_com(12,13) = -dot(EwL_GeAz, MRotor_DrTr);    
C_com(12,14) = -dot(EwL_GeAz, MRotor_B1F1);   
C_com(12,15) = -dot(EwL_GeAz, MRotor_B1E1);   
C_com(12,16) = -dot(EwL_GeAz, MRotor_B1F2);  
C_com(12,17) = -dot(EwL_GeAz, MRotor_B2F1);  
C_com(12,18) = -dot(EwL_GeAz, MRotor_B2E1);   
C_com(12,19) = -dot(EwL_GeAz, MRotor_B2F2);   
C_com(12,20) = -dot(EwL_GeAz, MRotor_B3F1);   
C_com(12,21) = -dot(EwL_GeAz, MRotor_B3E1);   
C_com(12,22) = -dot(EwL_GeAz, MRotor_B3F2);  

% r = 13  
C_com(13,13) = -dot(EwL_DrTr, MRotor_DrTr);    
C_com(13,14) = -dot(EwL_DrTr, MRotor_B1F1);   
C_com(13,15) = -dot(EwL_DrTr, MRotor_B1E1);   
C_com(13,16) = -dot(EwL_DrTr, MRotor_B1F2);  
C_com(13,17) = -dot(EwL_DrTr, MRotor_B2F1);  
C_com(13,18) = -dot(EwL_DrTr, MRotor_B2E1);   
C_com(13,19) = -dot(EwL_DrTr, MRotor_B2F2);   
C_com(13,20) = -dot(EwL_DrTr, MRotor_B3F1);   
C_com(13,21) = -dot(EwL_DrTr, MRotor_B3E1);   
C_com(13,22) = -dot(EwL_DrTr, MRotor_B3F2);  

C_com  = triu(C_com)+triu(C_com,1)';

%% Additional C terms
% Blade 1
C_B1(14,14) = trapz(BldSec, mB1.*dot(EvS1_B1F1,EvS1_B1F1,2), 3);
C_B1(14,15) = trapz(BldSec, mB1.*dot(EvS1_B1F1,EvS1_B1E1,2), 3);
C_B1(14,16) = trapz(BldSec, mB1.*dot(EvS1_B1F1,EvS1_B1F2,2), 3);

C_B1(15,15) = trapz(BldSec, mB1.*dot(EvS1_B1E1,EvS1_B1E1,2), 3);
C_B1(15,16) = trapz(BldSec, mB1.*dot(EvS1_B1E1,EvS1_B1F2,2), 3);

C_B1(16,16) = trapz(BldSec, mB1.*dot(EvS1_B1F2,EvS1_B1F2,2), 3);

C_B1  = triu(C_B1)+triu(C_B1,1)';
% end of Blade 1

% #7. Blade 2
C_B2(17,17) = trapz(BldSec, mB2.*dot(EvS2_B2F1,EvS2_B2F1,2), 3);
C_B2(17,18) = trapz(BldSec, mB2.*dot(EvS2_B2F1,EvS2_B2E1,2), 3);
C_B2(17,19) = trapz(BldSec, mB2.*dot(EvS2_B2F1,EvS2_B2F2,2), 3);

C_B2(18,18) = trapz(BldSec, mB2.*dot(EvS2_B2E1,EvS2_B2E1,2), 3);
C_B2(18,19) = trapz(BldSec, mB2.*dot(EvS2_B2E1,EvS2_B2F2,2), 3);

C_B2(19,19) = trapz(BldSec, mB2.*dot(EvS2_B2F2,EvS2_B2F2,2), 3);

C_B2  = triu(C_B2)+triu(C_B2,1)';
% end of Blade 2

% Blade 3
C_B3(20,20) = trapz(BldSec, mB3.*dot(EvS3_B3F1,EvS3_B3F1,2), 3);
C_B3(20,21) = trapz(BldSec, mB3.*dot(EvS3_B3F1,EvS3_B3E1,2), 3);
C_B3(20,22) = trapz(BldSec, mB3.*dot(EvS3_B3F1,EvS3_B3F2,2), 3);

C_B3(21,21) = trapz(BldSec, mB3.*dot(EvS3_B3E1,EvS3_B3E1,2), 3);
C_B3(21,22) = trapz(BldSec, mB3.*dot(EvS3_B3E1,EvS3_B3F2,2), 3);

C_B3(22,22) = trapz(BldSec, mB3.*dot(EvS3_B3F2,EvS3_B3F2,2), 3);

C_B3  = triu(C_B3)+triu(C_B3,1)';
% end of Blade 3

% Tower
C_T(7,7)   = trapz(TwrSec, mT.*dot(EvF_TFA1,EvF_TFA1,2), 3) + YawBrMass*dot(EvO_TFA1,EvO_TFA1); 
C_T(7,8)   = trapz(TwrSec, mT.*dot(EvF_TFA1,EvF_TSS1,2), 3) + YawBrMass*dot(EvO_TFA1,EvO_TSS1);
C_T(7,9)   = trapz(TwrSec, mT.*dot(EvF_TFA1,EvF_TFA2,2), 3) + YawBrMass*dot(EvO_TFA1,EvO_TFA2);
C_T(7,10)  = trapz(TwrSec, mT.*dot(EvF_TFA1,EvF_TSS2,2), 3) + YawBrMass*dot(EvO_TFA1,EvO_TSS2);

C_T(8,8)   = trapz(TwrSec, mT.*dot(EvF_TSS1,EvF_TSS1,2), 3) + YawBrMass*dot(EvO_TSS1,EvO_TSS1); 
C_T(8,9)   = trapz(TwrSec, mT.*dot(EvF_TSS1,EvF_TFA2,2), 3) + YawBrMass*dot(EvO_TSS1,EvO_TFA2);
C_T(8,10)  = trapz(TwrSec, mT.*dot(EvF_TSS1,EvF_TSS2,2), 3) + YawBrMass*dot(EvO_TSS1,EvO_TSS2);

C_T(9,9)   = trapz(TwrSec, mT.*dot(EvF_TFA2,EvF_TFA2,2), 3) + YawBrMass*dot(EvO_TFA2,EvO_TFA2);
C_T(9,10)  = trapz(TwrSec, mT.*dot(EvF_TFA2,EvF_TSS2,2), 3) + YawBrMass*dot(EvO_TFA2,EvO_TSS2);

C_T(10,10) = trapz(TwrSec, mT.*dot(EvF_TSS2,EvF_TSS2,2), 3) + YawBrMass*dot(EvO_TSS2,EvO_TSS2);

C_T  = triu(C_T)+triu(C_T,1)';
% end of Tower

%% Common f terms
f_com(1)  = dot(EvZ_Sg,   FAll_t);
f_com(2)  = dot(EvZ_Sw,   FAll_t);
f_com(3)  = dot(EvZ_Hv,   FAll_t);
f_com(4)  = dot(EwX_R,    MAll_t);
f_com(5)  = dot(EwX_P,    MAll_t);
f_com(6)  = dot(EwX_Y,    MAll_t);
f_com(7)  = dot(EvO_TFA1, FO_t) + dot(EwB_TFA1, MO_t);
f_com(8)  = dot(EvO_TSS1, FO_t) + dot(EwB_TSS1, MO_t);
f_com(9)  = dot(EvO_TFA2, FO_t) + dot(EwB_TFA2, MO_t);
f_com(10) = dot(EvO_TSS2, FO_t) + dot(EwB_TSS2, MO_t);
f_com(11) = dot(EwN_yaw,  MO_t);
f_com(12) = dot(EwL_GeAz, MRotor_t);
f_com(13) = dot(EwL_DrTr, MRotor_t);

%% Additional f terms
% Blade 1

f_B1(14) = trapz(BldSec, dot(EvS1_B1F1, F_AeroB1 - coprod(mB1,(g*repmat(Z(2,:),[1,1,nb]) + dEvS1)), 2) + dot(EwM1_B1F1, M_AeroB1, 2), 3);
f_B1(15) = trapz(BldSec, dot(EvS1_B1E1, F_AeroB1 - coprod(mB1,(g*repmat(Z(2,:),[1,1,nb]) + dEvS1)), 2) + dot(EwM1_B1E1, M_AeroB1, 2), 3);
f_B1(16) = trapz(BldSec, dot(EvS1_B1F2, F_AeroB1 - coprod(mB1,(g*repmat(Z(2,:),[1,1,nb]) + dEvS1)), 2) + dot(EwM1_B1F2, M_AeroB1, 2), 3);

f_ElasticB1(14) = -k11_B1F*q_B1F1 - k12_B1F*q_B1F2;
f_ElasticB1(15) = -k11_B1E*q_B1E1;
f_ElasticB1(16) = -k21_B1F*q_B1F1 - k22_B1F*q_B1F2;

f_DampingB1(14) = -zeta1_B1F*k11_B1F/pi/f1_B1F*qd_B1F1 - zeta2_B1F*k12_B1F/pi/f2_B1F*qd_B1F2;
f_DampingB1(15) = -zeta1_B1E*k11_B1E/pi/f1_B1E*qd_B1E1;
f_DampingB1(16) = -zeta1_B1F*k21_B1F/pi/f1_B1F*qd_B1F1 - zeta2_B1F*k22_B1F/pi/f2_B1F*qd_B1F2;

% Blade 2
f_B2(17) = trapz(BldSec, dot(EvS2_B2F1, F_AeroB2 - coprod(mB2,(g*repmat(Z(2,:),[1,1,nb]) + dEvS2)), 2) + dot(EwM2_B2F1, M_AeroB2, 2), 3);
f_B2(18) = trapz(BldSec, dot(EvS2_B2E1, F_AeroB2 - coprod(mB2,(g*repmat(Z(2,:),[1,1,nb]) + dEvS2)), 2) + dot(EwM2_B2E1, M_AeroB2, 2), 3);
f_B2(19) = trapz(BldSec, dot(EvS2_B2F2, F_AeroB2 - coprod(mB2,(g*repmat(Z(2,:),[1,1,nb]) + dEvS2)), 2) + dot(EwM2_B2F2, M_AeroB2, 2), 3);

f_ElasticB2(17) = -k11_B2F*q_B2F1 - k12_B2F*q_B2F2;
f_ElasticB2(18) = -k11_B2E*q_B2E1;
f_ElasticB2(19) = -k21_B2F*q_B2F1 - k22_B2F*q_B2F2;

f_DampingB2(17) = -zeta1_B2F*k11_B2F/pi/f1_B2F*qd_B2F1 - zeta2_B2F*k12_B2F/pi/f2_B2F*qd_B2F2;
f_DampingB2(18) = -zeta1_B2E*k11_B2E/pi/f1_B2E*qd_B2E1;
f_DampingB2(19) = -zeta1_B2F*k21_B2F/pi/f1_B2F*qd_B2F1 - zeta2_B2F*k22_B2F/pi/f2_B2F*qd_B2F2;

% Blade 3
f_B3(20) = trapz(BldSec, dot(EvS3_B3F1, F_AeroB3 - coprod(mB3,(g*repmat(Z(2,:),[1,1,nb]) + dEvS3)), 2) + dot(EwM3_B3F1, M_AeroB3, 2), 3);
f_B3(21) = trapz(BldSec, dot(EvS3_B3E1, F_AeroB3 - coprod(mB3,(g*repmat(Z(2,:),[1,1,nb]) + dEvS3)), 2) + dot(EwM3_B3E1, M_AeroB3, 2), 3);
f_B3(22) = trapz(BldSec, dot(EvS3_B3F2, F_AeroB3 - coprod(mB3,(g*repmat(Z(2,:),[1,1,nb]) + dEvS3)), 2) + dot(EwM3_B3F2, M_AeroB3, 2), 3);

f_ElasticB3(20) = -k11_B3F*q_B3F1 - k12_B3F*q_B3F2;
f_ElasticB3(21) = -k11_B3E*q_B3E1;
f_ElasticB3(22) = -k21_B3F*q_B3F1 - k22_B3F*q_B3F2;

f_DampingB3(20) = -zeta1_B3F*k11_B3F/pi/f1_B3F*qd_B3F1 - zeta2_B3F*k12_B3F/pi/f2_B3F*qd_B3F2;
f_DampingB3(21) = -zeta1_B3E*k11_B3E/pi/f1_B3E*qd_B3E1;
f_DampingB3(22) = -zeta1_B3F*k21_B3F/pi/f1_B3F*qd_B3F1 - zeta2_B3F*k22_B3F/pi/f2_B3F*qd_B3F2;

% Tower
f_T(7)  = trapz(TwrSec, dot(EvF_TFA1, F_AeroT - coprod(mT,(g*repmat(Z(2,:),[1,1,nt]) + dEvF)), 2), 3) - YawBrMass*dot(EvO_TFA1, g*Z(2,:) + dEvO);
f_T(8)  = trapz(TwrSec, dot(EvF_TSS1, F_AeroT - coprod(mT,(g*repmat(Z(2,:),[1,1,nt]) + dEvF)), 2), 3) - YawBrMass*dot(EvO_TSS1, g*Z(2,:) + dEvO);
f_T(9)  = trapz(TwrSec, dot(EvF_TFA2, F_AeroT - coprod(mT,(g*repmat(Z(2,:),[1,1,nt]) + dEvF)), 2), 3) - YawBrMass*dot(EvO_TFA2, g*Z(2,:) + dEvO);
f_T(10) = trapz(TwrSec, dot(EvF_TSS2, F_AeroT - coprod(mT,(g*repmat(Z(2,:),[1,1,nt]) + dEvF)), 2), 3) - YawBrMass*dot(EvO_TSS2, g*Z(2,:) + dEvO);

f_ElasticT(7)  = -k11_TFA*q_TFA1 - k12_TFA*q_TFA2;
f_ElasticT(8)  = -k11_TSS*q_TSS1 - k12_TSS*q_TSS2;
f_ElasticT(9)  = -k21_TFA*q_TFA1 - k22_TFA*q_TFA2;
f_ElasticT(10) = -k21_TSS*q_TSS1 - k22_TSS*q_TSS2;

f_DampingT(7)  = -zeta1_TFA*k11_TFA/pi/f1_TFA*qd_TFA1 - zeta2_TFA*k12_TFA/pi/f2_TFA*qd_TFA2;
f_DampingT(8)  = -zeta1_TSS*k11_TSS/pi/f1_TSS*qd_TSS1 - zeta2_TSS*k12_TSS/pi/f2_TSS*qd_TSS2;
f_DampingT(9)  = -zeta1_TFA*k21_TFA/pi/f1_TFA*qd_TFA1 - zeta2_TFA*k22_TFA/pi/f2_TFA*qd_TFA2;
f_DampingT(10) = -zeta1_TSS*k21_TSS/pi/f1_TSS*qd_TSS1 - zeta2_TSS*k22_TSS/pi/f2_TSS*qd_TSS2;


% Yaw spring and damper
f_Yaw(11) = -YawSpr*(q_yaw - YawNeut) - YawDamp*qd_yaw;

% Generator & Drive-train
f_G(12)            = - dot(EwG_GeAz, (I_G*dEwN')');
T_Brake            = 0.0;
f_Gen(12)          = -GBRatio*GenTrq;
f_Brake(12)        = -GBRatio*T_Brake;
f_ElasticDrive(13) = -DTTorSpr*q_DrTr;
f_DampDrive(13)    = -DTTorDmp*qd_DrTr;

%% Update additional output vector
Controls(5:3*nb+4) = reshape(phi, [1 3*nb]);

%% Assembly
IM_nom  = C_com + C_B1 + C_B2 + C_B3 + C_T;
IM_nom(1:6,1:6) = IM_nom(1:6,1:6) + Platform.AM;

f_nom = f_com + f_B1 + f_B2 + f_B3 + f_ElasticB1 + f_DampingB1 + f_ElasticB2 + f_DampingB2 +f_ElasticB3 + f_DampingB3 + f_T + f_ElasticT + f_DampingT + f_Yaw  + f_Gen + f_Brake ...
    + f_ElasticDrive + f_DampDrive + f_G;
f_nom(1:6) = f_nom(1:6) + f_AddDamp;

f_nom(6)   = f_nom(6) - 98340000*(q_Y - (-1.2*10^-3));


end

