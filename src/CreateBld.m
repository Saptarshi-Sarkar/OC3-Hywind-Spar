function Bld = CreateBld(ElastoDyn,Geometry,Blade)
% Blade shape functions
BldLen = (ElastoDyn.TipRad-ElastoDyn.HubRad);
Bld.BldSec = Geometry.Blade(:,1);
Bld.Twist   = Geometry.Blade(:,5);
Bld.nb   = length(Bld.BldSec);

for i = 1:length(Blade.FlapMode1)
    Blade.FlapMode1(i) = Blade.FlapMode1(i)/BldLen^(length(Blade.FlapMode1)-i);
    Blade.FlapMode2(i) = Blade.FlapMode2(i)/BldLen^(length(Blade.FlapMode2)-i);
    Blade.EdgeMode1(i) = Blade.EdgeMode1(i)/BldLen^(length(Blade.EdgeMode1)-i);
end

[O1_B1, O2_B1, O3_B1, W1_B1, W2_B1, W3_B1, dO1_B1, dO2_B1, dO3_B1, dW1_B1, dW2_B1,...
    dW3_B1, s11_B1, s22_B1, s33_B1, s12_B1, s13_B1, s23_B1 ] ...
    = BladeModeShapes( Blade.FlapMode1, Blade.FlapMode2, Blade.EdgeMode1, Bld.Twist, Bld.BldSec );

Bld.O1_B1 = reshape(O1_B1,1,1,Bld.nb);
Bld.O2_B1 = reshape(O2_B1,1,1,Bld.nb);
Bld.O3_B1 = reshape(O3_B1,1,1,Bld.nb);

Bld.W1_B1 = reshape(W1_B1,1,1,Bld.nb);
Bld.W2_B1 = reshape(W2_B1,1,1,Bld.nb);
Bld.W3_B1 = reshape(W3_B1,1,1,Bld.nb);

Bld.s11_B1 = reshape(s11_B1,1,1,Bld.nb);
Bld.s22_B1 = reshape(s22_B1,1,1,Bld.nb);
Bld.s33_B1 = reshape(s33_B1,1,1,Bld.nb);

Bld.s12_B1 = reshape(s12_B1,1,1,Bld.nb);
Bld.s13_B1 = reshape(s13_B1,1,1,Bld.nb);
Bld.s23_B1 = reshape(s23_B1,1,1,Bld.nb);

ddBlFlapMode1 = polyval(polyder(polyder(Blade.FlapMode1)),Blade.Data(:,1)*BldLen);
ddBlFlapMode2 = polyval(polyder(polyder(Blade.FlapMode2)),Blade.Data(:,1)*BldLen);
ddBlEdgeMode1 = polyval(polyder(polyder(Blade.EdgeMode1)),Blade.Data(:,1)*BldLen);

Blade.Data(:,4) = Blade.Data(:,4)*1.0456;

Bld.mB1     = reshape(interp1(Blade.Data(:,1)*BldLen, Blade.Data(:,4), Bld.BldSec),1,1,Bld.nb); 

Bld.k11_B1F = simps(Blade.Data(:,1)*BldLen, Blade.Data(:,5).*ddBlFlapMode1.*ddBlFlapMode1);
Bld.k12_B1F = simps(Blade.Data(:,1)*BldLen, Blade.Data(:,5).*ddBlFlapMode1.*ddBlFlapMode2);
Bld.k22_B1F = simps(Blade.Data(:,1)*BldLen, Blade.Data(:,5).*ddBlFlapMode2.*ddBlFlapMode2);

Bld.k11_B1E = simps(Blade.Data(:,1)*BldLen, Blade.Data(:,6).*ddBlEdgeMode1.*ddBlEdgeMode1);

m11_B1F = simps(Blade.Data(:,1)*BldLen, Blade.Data(:,4).*polyval(conv(Blade.FlapMode1,Blade.FlapMode1), Blade.Data(:,1)*BldLen));
m22_B1F = simps(Blade.Data(:,1)*BldLen, Blade.Data(:,4).*polyval(conv(Blade.FlapMode2,Blade.FlapMode2), Blade.Data(:,1)*BldLen));
m11_B1E = simps(Blade.Data(:,1)*BldLen, Blade.Data(:,4).*polyval(conv(Blade.EdgeMode1,Blade.EdgeMode1), Blade.Data(:,1)*BldLen));

Bld.f1_B1F = sqrt(Bld.k11_B1F/m11_B1F)/2/pi;
Bld.f2_B1F = sqrt(Bld.k22_B1F/m22_B1F)/2/pi;
Bld.f1_B1E = sqrt(Bld.k11_B1E/m11_B1E)/2/pi;

Bld.zeta1_B1F = Blade.DampBF1/100;
Bld.zeta2_B1F = Blade.DampBF2/100;
Bld.zeta1_B1E = Blade.DampBE1/100;

% These quantities are of load calculation
% Load.BlSpn   = Geometry.Blade(:,1);
Bld.Twist      = Geometry.Blade(:,5);
Bld.Cord       = Geometry.Blade(:,6);
Bld.FoilNo     = Geometry.Blade(:,7);
Bld.PitchAxis  = interp1(Blade.Data(:,1)*61.5,Blade.Data(:,2), Geometry.Blade(:,1));
TempDist       = (0.25 - Bld.PitchAxis).*Bld.Cord;
TempDistJ1     = TempDist.*sin(Bld.Twist);
TempDistJ2     = TempDist.*cos(Bld.Twist);
Bld.AeroCentJ1 = TempDistJ1.*cos(Bld.Twist) - TempDistJ2.*sin(Bld.Twist);
Bld.AeroCentJ2 = TempDistJ1.*sin(Bld.Twist) + TempDistJ2.*cos(Bld.Twist);
Bld.dO1_B1     = reshape(dO1_B1,1,1,Bld.nb);  
Bld.dO2_B1     = reshape(dO2_B1,1,1,Bld.nb); 
Bld.dO3_B1     = reshape(dO3_B1,1,1,Bld.nb); 

Bld.dW1_B1     = reshape(dW1_B1,1,1,Bld.nb); 
Bld.dW2_B1     = reshape(dW2_B1,1,1,Bld.nb); 
Bld.dW3_B1     = reshape(dW3_B1,1,1,Bld.nb);        %hold on; plot(Load.BlSpn,reshape(Load.dO1_B1,Load.nb_Load,1))

end

