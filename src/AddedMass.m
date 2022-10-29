function [ AM ] = AddedMass( Platform )
global FLUIDDENSITY;
CA = 0.969954;
ndz         = length(Platform.StripDepths);
StripDia    = Platform.StripDia;
StripDepths = Platform.StripDepths;
Z1  = repmat([1 0 0],1,1,ndz);
Z3  = repmat([0 0 1],1,1,ndz);
rZZ = reshape(-StripDepths',1,1,ndz).*[0 1 0];
AM     = zeros(6,6);

%% Added Mass Coeffects

AM(1,1) = CA*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDia.^2);
AM(1,5) = -CA*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDia.^2.* reshape(dot(Z1, cross(Z3, rZZ, 2), 2),1, ndz));

AM(2,2) = CA*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDia.^2);
AM(2,4) = -CA*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDia.^2.* reshape(dot(Z3, cross(Z1, rZZ, 2), 2), 1, ndz));

AM(3,3) = FLUIDDENSITY*pi*StripDia(end)^3/12;

AM(4,2) = CA*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDia.^2.*StripDepths);
AM(4,4) = -CA*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDia.^2.*StripDepths.*reshape(dot(Z3, cross(Z1, rZZ, 2), 2), 1, ndz));

AM(5,1) = -CA*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDia.^2.*StripDepths);
AM(5,5) = CA*FLUIDDENSITY*pi/4*trapz(StripDepths, StripDia.^2.*StripDepths.*reshape(dot(Z1, cross(Z3, rZZ, 2), 2), 1, ndz));


end

