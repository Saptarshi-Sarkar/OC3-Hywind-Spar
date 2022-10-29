function [IM, f, Controls] = SystemMatrices(t, q, Controls, DOFs, ElastoDyn, Airfoils, Twr, Bld, Platform, Wind, WindNom, wave, mooring_load_ptr, Servo)
%#Codegen
ep = 1e-15;
OnePlusEps = 1+ep;
persistent LastTime
if isempty(LastTime)
    LastTime = 0;
    fprintf('Persistent variables in System Matrices initialized.\n');
end

%% Update state
% Dofs
Comq = zeros(DOFs.Avail*2,1);
IndexActive = [DOFs.Active, DOFs.Active+DOFs.Avail];
IndexNominal = [DOFs.ActNominal, DOFs.ActNominal+DOFs.Avail];

Comq(IndexActive) = q';
q_Nom  = zeros(44,1);
q_Nom([DOFs.ActNominal, DOFs.ActNominal+22]) = Comq(IndexNominal);
q_Sg   = Comq(1);
q_Sw   = Comq(2);
q_Hv   = Comq(3);
q_R    = Comq(4);
q_P    = Comq(5);
q_Y    = Comq(6);
% q_TFA1 = Comq(7);
% q_TSS1 = Comq(8);
% q_TFA2 = Comq(9);
% q_TSS2 = Comq(10);
% q_yaw  = Comq(11);
% q_GeAz = Comq(12);
% q_DrTr = Comq(13);
% q_B1F1 = Comq(14);
% q_B1E1 = Comq(15);
% q_B1F2 = Comq(16);
% q_B2F1 = Comq(17);
% q_B2E1 = Comq(18);
% q_B2F2 = Comq(19);
% q_B3F1 = Comq(20);
% q_B3E1 = Comq(21);
% q_B3F2 = Comq(22);
 
qd_Sg   = Comq(24);
qd_Sw   = Comq(25);
qd_Hv   = Comq(26);
qd_R    = Comq(27);
qd_P    = Comq(28);
qd_Y    = Comq(29);
% qd_TFA1 = Comq(30);
% qd_TSS1 = Comq(31);
% qd_TFA2 = Comq(32);
% qd_TSS2 = Comq(33);
% qd_yaw  = Comq(34);
qd_GeAz = Comq(35);

% qd_DrTr = Comq(36);
% qd_B1F1 = Comq(37);
% qd_B1E1 = Comq(38);
% qd_B1F2 = Comq(39);
% qd_B2F1 = Comq(40);
% qd_B2E1 = Comq(41);
% qd_B2F2 = Comq(42);
% qd_B3F1 = Comq(43);
% qd_B3E1 = Comq(44);
% qd_B3F2 = Comq(45);

% end of Update state

%% Call baseline controllers
BlPitch = Controls(2:4);
GenSpeed = ElastoDyn.GBRatio*(qd_GeAz);
if WindNom.PitchControl
    [GenTrq, BlPitch] = BaselineControllers(t, GenSpeed, BlPitch, Servo);
else
    [GenTrq, ~] = BaselineControllers(t, GenSpeed, BlPitch, Servo);
end

%% Update additional output vector
Controls(1:4) = [GenTrq, BlPitch];

%% Mooring loads
X  = [q_Sg, q_Sw, q_Hv, q_R, q_P, q_Y];
XD = [qd_Sg, qd_Sw, qd_Hv, qd_R, qd_P, qd_Y];

dt = t*OnePlusEps - LastTime; if dt < 0; disp(LastTime); disp(t); error('Time step negative'); end
if t == 0 || dt >= 1.25e-2
    if Platform.Mooring == 1
        calllib('MoorDyn', 'LinesCalc', X, XD, mooring_load_ptr, t, dt);
    elseif Platform.Mooring == 2
        calllib('MoorApiwin64','update', mooring_load_ptr, X', XD', t, dt, 0); 
    end
    LastTime = t;
end

mooring_load = mooring_load_ptr.value;

Z = eye(3); 
EwX = qd_R*Z(1,:) + qd_Y*Z(2,:) - qd_P*Z(3,:);
EvZ   = qd_Sg*Z(1,:) + qd_Hv*Z(2,:) - qd_Sw*Z(3,:);

f_Morison = Morisons(t, Z, EwX, EvZ, wave, Platform, Platform.WaveLoads);

%% Get the nominal system matrixes
WindNom.velocity = zeros(3,31,31);

qr = {t, [1 2 3], Wind.y, Wind.z};
WindNom.velocity(:) = Wind.Velocity(qr);

[IM_nom, f_nom, Controls] = NominalSystemMatrix_mex(q_Nom, Controls, ElastoDyn, Airfoils, Twr, Bld, Platform, WindNom, mooring_load, f_Morison);

%% Assemble Initial Inertia matrix and force vector
IM = zeros(DOFs.Avail,DOFs.Avail);
f  = zeros(DOFs.Avail,1);

IM(1:22,1:22) = IM_nom;
f(1:22,1)     = f_nom;

IM = IM(DOFs.Active,DOFs.Active);

f = f(DOFs.Active);

end
