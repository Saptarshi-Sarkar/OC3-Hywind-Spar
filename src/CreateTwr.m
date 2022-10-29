function Twr = CreateTwr(Tower,ElastoDyn)
% Tower mode shapes
TwrFlexHt            = ElastoDyn.TwrHt - ElastoDyn.TowerBsHt;
Twr.nt = Tower.NoOfPoints;
for i = 1:length(Tower.FAMode1)
    Tower.FAMode1(i) = Tower.FAMode1(i)/TwrFlexHt^(length(Tower.FAMode1)-i);
    Tower.FAMode2(i) = Tower.FAMode2(i)/TwrFlexHt^(length(Tower.FAMode2)-i);
    Tower.SSMode1(i) = Tower.SSMode1(i)/TwrFlexHt^(length(Tower.SSMode1)-i);
    Tower.SSMode2(i) = Tower.SSMode2(i)/TwrFlexHt^(length(Tower.SSMode2)-i);
end

Twr.TwrSec = Tower.Data(:,1)*TwrFlexHt;
Twr.O1_TFA = reshape(polyval(Tower.FAMode1,Twr.TwrSec),1,1,Twr.nt);             %plot(Twr.TwrSec,reshape(Twr.O1_TFA,Twr.nt,1))
Twr.O2_TFA = reshape(polyval(Tower.FAMode2,Twr.TwrSec),1,1,Twr.nt);             %hold on; plot(Twr.TwrSec,reshape(Twr.O2_TFA,Twr.nt,1))
Twr.O1_TSS = reshape(polyval(Tower.SSMode1,Twr.TwrSec),1,1,Twr.nt);             %plot(TwrSec,reshape(O1_TSS,nt,1))
Twr.O2_TSS = reshape(polyval(Tower.SSMode2,Twr.TwrSec),1,1,Twr.nt);             %plot(TwrSec,reshape(O2_TSS,nt,1))

dFAMode1 = polyder(Tower.FAMode1);
dFAMode2 = polyder(Tower.FAMode2);
dSSMode1 = polyder(Tower.SSMode1);
dSSMode2 = polyder(Tower.SSMode2);
ddFAMode1 = polyder(dFAMode1);
ddFAMode2 = polyder(dFAMode2);
ddSSMode1 = polyder(dSSMode1);
ddSSMode2 = polyder(dSSMode2);

Twr.s11_TFA = reshape(polyval(polyint(conv(dFAMode1,dFAMode1)), Twr.TwrSec),1,1,Twr.nt);
Twr.s12_TFA = reshape(polyval(polyint(conv(dFAMode1,dFAMode2)), Twr.TwrSec),1,1,Twr.nt);  
Twr.s22_TFA = reshape(polyval(polyint(conv(dFAMode2,dFAMode2)), Twr.TwrSec),1,1,Twr.nt);

Twr.S11_TFA = polyval(polyint(conv(dFAMode1,dFAMode1)), TwrFlexHt);             
Twr.S12_TFA = polyval(polyint(conv(dFAMode1,dFAMode2)), TwrFlexHt);
Twr.S22_TFA = polyval(polyint(conv(dFAMode2,dFAMode2)), TwrFlexHt);

Twr.s11_TSS = reshape(polyval(polyint(conv(dSSMode1,dSSMode1)), Twr.TwrSec),1,1,Twr.nt);
Twr.s12_TSS = reshape(polyval(polyint(conv(dSSMode1,dSSMode2)), Twr.TwrSec),1,1,Twr.nt);
Twr.s22_TSS = reshape(polyval(polyint(conv(dSSMode2,dSSMode2)), Twr.TwrSec),1,1,Twr.nt);

Twr.S11_TSS = polyval(polyint(conv(dSSMode1,dSSMode1)), TwrFlexHt);             
Twr.S12_TSS = polyval(polyint(conv(dSSMode1,dSSMode2)), TwrFlexHt);
Twr.S22_TSS = polyval(polyint(conv(dSSMode2,dSSMode2)), TwrFlexHt);

Twr.dO1_TFA = polyval(polyder(Tower.FAMode1),TwrFlexHt);
Twr.dO1_TSS = polyval(polyder(Tower.SSMode1),TwrFlexHt);
Twr.dO2_TFA = polyval(polyder(Tower.FAMode2),TwrFlexHt);
Twr.dO2_TSS = polyval(polyder(Tower.SSMode2),TwrFlexHt);

Twr.mT      = Tower.Data(:,2); TwrMass = simps(Twr.TwrSec,Twr.mT);

Twr.k11_TFA = simps(Twr.TwrSec,Tower.Data(:,3).*polyval(conv(ddFAMode1,ddFAMode1),Twr.TwrSec));           % when lower limit is zero, taking limits and evaluating 
Twr.k12_TFA = simps(Twr.TwrSec,Tower.Data(:,3).*polyval(conv(ddFAMode1,ddFAMode2),Twr.TwrSec));           % the value at final point is same. Here two methods are used
Twr.k11_TSS = simps(Twr.TwrSec,Tower.Data(:,4).*polyval(conv(ddSSMode1,ddSSMode1),Twr.TwrSec));           % to make the point
Twr.k12_TSS = simps(Twr.TwrSec,Tower.Data(:,4).*polyval(conv(ddSSMode1,ddSSMode2),Twr.TwrSec));
Twr.k21_TFA = simps(Twr.TwrSec,Tower.Data(:,3).*polyval(conv(ddFAMode2,ddFAMode1),Twr.TwrSec));
Twr.k22_TFA = simps(Twr.TwrSec,Tower.Data(:,3).*polyval(conv(ddFAMode2,ddFAMode2),Twr.TwrSec));
Twr.k21_TSS = simps(Twr.TwrSec,Tower.Data(:,4).*polyval(conv(ddSSMode2,ddSSMode1),Twr.TwrSec));
Twr.k22_TSS = simps(Twr.TwrSec,Tower.Data(:,4).*polyval(conv(ddSSMode2,ddSSMode2),Twr.TwrSec));

m11_TFA = simps(Twr.TwrSec,Twr.mT.*polyval(conv(Tower.FAMode1,Tower.FAMode1),Twr.TwrSec));
m22_TFA = simps(Twr.TwrSec,Twr.mT.*polyval(conv(Tower.FAMode2,Tower.FAMode2),Twr.TwrSec));
m11_TSS = simps(Twr.TwrSec,Twr.mT.*polyval(conv(Tower.SSMode1,Tower.SSMode1),Twr.TwrSec));
m22_TSS = simps(Twr.TwrSec,Twr.mT.*polyval(conv(Tower.SSMode2,Tower.SSMode2),Twr.TwrSec));

m_dummy_TFA1 = m11_TFA + ElastoDyn.HubMass + ElastoDyn.NacMass + 17740*3 ;
m_dummy_TFA2 = m22_TFA + ElastoDyn.HubMass + ElastoDyn.NacMass + 17740*3 ;
m_dummy_TSS1 = m11_TSS + ElastoDyn.HubMass + ElastoDyn.NacMass + 17740*3 ;

Twr.f1_TFA = sqrt(Twr.k11_TFA/m11_TFA)/2/pi;
Twr.f2_TFA = sqrt(Twr.k22_TFA/m22_TFA)/2/pi;
Twr.f1_TSS = sqrt(Twr.k11_TSS/m11_TSS)/2/pi;
Twr.f2_TSS = sqrt(Twr.k22_TSS/m22_TSS)/2/pi;

Twr.zeta1_TFA = Tower.Damp_TFA1/100;
Twr.zeta2_TFA = Tower.Damp_TFA2/100;
Twr.zeta1_TSS = Tower.Damp_TSS1/100;
Twr.zeta2_TSS = Tower.Damp_TSS2/100;

Twr.mT      = reshape(Tower.Data(:,2),1,1,Twr.nt);
end

