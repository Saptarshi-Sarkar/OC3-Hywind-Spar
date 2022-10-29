function [ComGenTrq, ComBlPitch] = BaselineControllers( t, GenSpeed, BlPitch, Servo )
%#Codegen
persistent GenSpeedF IntSpdErr LastGenTrq LastTime LastTimePC LastTimeVS PitCom

OnePlusEps = 1+eps;
if isempty(GenSpeedF)
    if isempty(IntSpdErr)~=1 || isempty(LastGenTrq)~=1 || isempty(LastTime)~=1 || isempty(LastTimePC)~=1 || ...
            isempty(LastTimeVS)~=1 || isempty(PitCom)~=1 
        error('All variables are not empty')
    end
    GenSpeedF = GenSpeed;
    PitCom = BlPitch;
    GK = 1/(1+PitCom(1)/Servo.PC_KK);
    IntSpdErr  = PitCom(1)/(GK*Servo.PC_KI);
    LastTime   = t;
    LastTimePC = t - Servo.PC_DT;
    LastTimeVS = t - Servo.VS_DT;
    LastGenTrq = Servo.VS_RtTq;
    fprintf('Persistent variables in Baseline Controllers initialized.\n');
end
PitRate = zeros(1,3);
%% Filter the HSS (generator) speed measurement
Alpha = exp((LastTime-t)*Servo.CornerFreq);
% Apply the filter
GenSpeedF = (1.0 - Alpha)*GenSpeed + Alpha*GenSpeedF;

%% Variable-speed torque control:
ElapTime = t - LastTimeVS;

if (t*OnePlusEps - LastTimeVS) >= Servo.VS_DT

    if GenSpeedF >= Servo.VS_RtGnSp || PitCom(1) >= Servo.VS_Rgn3MP  % Region 3
        GenTrq = Servo.VS_RtPwr/Servo.PC_RefSpd; 
    elseif GenSpeedF <= Servo.VS_CtInSp                                 % Region 1
        GenTrq = 0; 
    elseif GenSpeedF < Servo.VS_Rgn2Sp                                  % Region 1.5
        GenTrq = Servo.VS_Slope15*(GenSpeedF - Servo.VS_CtInSp); 
    elseif GenSpeedF < Servo.VS_TrGnSp                                  % Region 2
        GenTrq = Servo.VS_Rgn2K*GenSpeedF*GenSpeedF;
    else                                                                % Region 2.5
        GenTrq = Servo.VS_Slope25*(GenSpeedF - Servo.VS_SySp);
    end

    % Saturate the command torque using maximum torque limit
    GenTrq = min(GenTrq, Servo.VS_MaxTq);

    % Saturate the command torque using the torque rate limit
    TrqRate = (GenTrq - LastGenTrq)/ElapTime;
    TrqRate = min(max(TrqRate, -Servo.VS_MaxRat), Servo.VS_MaxRat);
    GenTrq = LastGenTrq + TrqRate*ElapTime;
    
    % Reset the valuse of LastTimeVS and LastGenTrq to the current values
    LastTimeVS = t;
    LastGenTrq = GenTrq;
end
ComGenTrq = LastGenTrq;

%% Pitch Control
ElapTime = t - LastTimePC;

if (t*OnePlusEps - LastTimePC) >= Servo.PC_DT 
    % Compute the gain scheduling correction factor based on previously
    % commanded pitch angle for Blade 1
    GK = 1/(1+PitCom(1)/Servo.PC_KK);
    
    % Compute the curret speed error and its integral w.r.t. time; saturate
    % the integral term using the pitch angle limits
    
    SpdErr = GenSpeedF - Servo.PC_RefSpd; %disp(GenSpeedF); 
    IntSpdErr = IntSpdErr + SpdErr*ElapTime;
    IntSpdErr = min(max(IntSpdErr, Servo.PC_MinPit/(GK*Servo.PC_KI)), Servo.PC_MaxPit/(GK*Servo.PC_KI));
    
    % Compute the pitch commands associated with the proportional and
    % integral gains:
    
    PitComP = GK*Servo.PC_KP * SpdErr;
    PitComI = GK*Servo.PC_KI * IntSpdErr;
    
    % Superimpose the individual commands to get the total pitch command
    PitComT = PitComP + PitComI;
    PitComT = min(max(PitComT, Servo.PC_MinPit), Servo.PC_MaxPit);

    % Saturate the overall commanded pitch using the pitch rate limit
    for i = 1:3
        PitRate(i) = (PitComT - BlPitch(i))/ElapTime;
        PitRate(i) = min(max(PitRate(i), -Servo.PC_MaxRat), Servo.PC_MaxRat);
        PitCom(i)  = BlPitch(i) + PitRate(i)*ElapTime;       
        PitCom(i)  = min( max( PitCom(i), Servo.PC_MinPit ), Servo.PC_MaxPit );
    end
    
    % Reset the value of LastTimePC to current value
    LastTimePC = t;
end
    % Collective pitch control
    ComBlPitch = [PitCom(1), PitCom(1), PitCom(1)];
    LastTime = t;
end

