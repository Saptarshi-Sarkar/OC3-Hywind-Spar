function [ Servo ] = ServoDyn()
fid = fopen('./5MW_Baseline/NRELOffshrBsline5MW_Onshore_ServoDyn.dat','r');
RPM2RPS = pi/30;

% Simple varible power torque control data

Info = textscan(fid,'%f %s',4,'delimiter','\n','HeaderLines',27);
% Servo.VS_RtGnSp = Info{1,1}(1)*RPM2RPS;             
Servo.VS_RtGnSp = 121.6805;                       % For Baseline Controller
Servo.VS_RtTq = Info{1,1}(2); Servo.VS_RtTq = 43093.55;
Servo.VS_Rgn2K = Info{1,1}(3)/(RPM2RPS*RPM2RPS);
Servo.VS_SlPc = 10;                               % For Baseline Controller Use this value 
% Servo.VS_SlPc = Info{1,1}(4);
Servo.VS_SySp = Servo.VS_RtGnSp/(1+.01*Servo.VS_SlPc);
Servo.VS_CtInSp = 70.16224; % in rads/s
Servo.VS_DT = 0.00125;
Servo.VS_MaxRat = 15000;
Servo.VS_MaxTq  = 47402.91;                         % 10% more than Rated Value
Servo.VS_Rgn2Sp = 91.21091; % in rads/s
Servo.VS_Rgn3MP = 0.01745329;
Servo.VS_RtPwr  = 5296610.0;

Servo.PC_RefSpd = 122.9096;

% ============== For Baseline Controller ==============
Servo.VS_SySp = Servo.VS_RtGnSp/(1+.01*Servo.VS_SlPc);
Servo.VS_Slope15 = (Servo.VS_Rgn2K*Servo.VS_Rgn2Sp*Servo.VS_Rgn2Sp)/(Servo.VS_Rgn2Sp - Servo.VS_CtInSp);
Servo.VS_Slope25 = (Servo.VS_RtPwr/Servo.PC_RefSpd)/ ( Servo.VS_RtGnSp - Servo.VS_SySp);

if Servo.VS_Rgn2K == 0
    Servo.VS_TrGnSp = Servo.VS_SySp;
else
    Servo.VS_TrGnSp = ( Servo.VS_Slope25 - sqrt( Servo.VS_Slope25*( Servo.VS_Slope25 - 4.0*Servo.VS_Rgn2K*Servo.VS_SySp ) ) )/( 2.0*Servo.VS_Rgn2K );
end
% ================= X ====== X =====================

% Servo.VS_Slope15 = (Servo.VS_Rgn2K*Servo.VS_Rgn2Sp*Servo.VS_Rgn2Sp)/(Servo.VS_Rgn2Sp - Servo.VS_CtInSp);
% if Servo.VS_SlPc < sqrt(eps(Servo.VS_SlPc)) % Then, We don't have a 2.5 so we'll use  VS_TrGnSp = VS_RtGnSp
%     disp('We donot have a 2.5 region')
%     Servo.VS_Slope25 = 9999.9;
%     Servo.VS_TrGnSp = Servo.VS_RtGnSp;
% else
%     Servo.VS_Slope25 = Servo.VS_RtTq /( Servo.VS_RtGnSp - Servo.VS_SySp ); % Torque/speed slope of region 2 1/2 induction generator.
%     if  abs(Servo.VS_Rgn2K) < eps(Servo.VS_SlPc) % TRUE. if the Region 2 torque is flat, and thus, the denominator in the ELSE condition is zero
%         disp('Region 2 torque is flat')
%         Servo.VS_TrGnSp = Servo.VS_SySp;
%     else
%         Servo.VS_TrGnSp = (Servo.VS_Slope25 - sqrt(Servo.VS_Slope25*(Servo.VS_Slope25 - 4*Servo.VS_Rgn2K*Servo.VS_SySp)))/(2*Servo.VS_Rgn2K);  % Transitional generator speed between regions 2 and 2 1/2.
%     end
% end

Servo.PC_DT = 0.00125;
Servo.PC_KI = 0.0008965149;
Servo.PC_KK = 0.1099965;
Servo.PC_KP = 0.006275604;
Servo.PC_MaxPit = 1.570796;
Servo.PC_MaxRat = 0.1396263;
Servo.PC_MinPit = 0;


Servo.CornerFreq = 2*pi*.25;

end

