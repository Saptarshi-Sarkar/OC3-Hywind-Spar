function [f_Morison] = Morisons(t, Z, EWX, EVZ, wave, Platform, Switch)
% Initializations
global FLUIDDENSITY;
CA = 0.969954;
CD = 0.6;
ndz         = length(Platform.StripDepths);
StripDia    = Platform.StripDia;
StripDepths = Platform.StripDepths;
Z1  = repmat(Z(1,:),1,1,ndz);
Z3  = repmat(Z(3,:),1,1,ndz);
rZZ = reshape(-StripDepths',1,1,ndz).*Z(2,:);
EwX = repmat(EWX,1,1,ndz);
EvZ = repmat(EVZ,1,1,ndz);
f_Morison = zeros(1,6);

%% Hydrodynamic forces and moments
EvZ1 = reshape(dot(Z1, EvZ + cross(EwX, rZZ, 2) , 2),1, ndz); %disp(EvZ1);
EvZ3 = reshape(dot(-Z3, EvZ + cross(EwX, rZZ, 2) , 2),1, ndz);

if Switch == 1
    wave.update_from_dataset(t);
    fluid_vel = wave.flow_velocity;
    fluid_acc = wave.flow_acceleration;
elseif Switch == 0
    fluid_vel = zeros(2, length(StripDepths));
    fluid_acc = zeros(2, length(StripDepths));
else
    error('Something is wrong with wave switch');
end

WaveAcc1 = fluid_acc(1,:);
WaveAcc3 = fluid_acc(2,:);
WaveVel1 = fluid_vel(1,:);
WaveVel3 = fluid_vel(2,:);

%% Forces and moments estimated per unit length and then integrated
% Forces
F_Surge = + (1+CA)*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDia.^2.*WaveAcc1) ...
          + CD*FLUIDDENSITY/2*trapz(StripDepths, StripDia.*(WaveVel1 - EvZ1).*abs(WaveVel1 - EvZ1));
           
F_Sway  = + (1+CA)*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDia.^2.*WaveAcc3)...
          + CD*FLUIDDENSITY/2*trapz(StripDepths, StripDia.*(WaveVel3 - EvZ3).*abs(WaveVel3 - EvZ3));
           
f_Morison(1) = F_Surge;
f_Morison(2) = F_Sway;

% Moments
M_Pitch = + (1+CA)*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDepths.*StripDia.^2.*WaveAcc1) ...
          + CD*FLUIDDENSITY/2*trapz(StripDepths, StripDepths.*StripDia.*(WaveVel1 - EvZ1).*abs(WaveVel1 - EvZ1));
          
M_Roll  = + (1+CA)*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDepths.*StripDia.^2.*WaveAcc3) ...
          + CD*FLUIDDENSITY/2*trapz(StripDepths, StripDepths.*StripDia.*(WaveVel3 - EvZ3).*abs(WaveVel3 - EvZ3));
          
f_Morison(4) = M_Roll;
f_Morison(5) = -M_Pitch;

end

