function phi = InitializeInflowAngle(ElastoDyn,Bld,Wind)
z   = zeros(Bld.nb, 3);
y   = zeros(Bld.nb, 3);
phi = zeros(Bld.nb, 3);
for iB = 1:3
    z(:,iB) = (ElastoDyn.HubRad + Bld.BldSec)*cosd((iB-1)*360/3) + ElastoDyn.TwrHt;
    y(:,iB) = (ElastoDyn.HubRad + Bld.BldSec)*sind((iB-1)*360/3);
end
Vy  = ElastoDyn.RotSpeed*(ElastoDyn.HubRad + Bld.BldSec);
for iB = 1:3
    for iNd = 1:Bld.nb
        Vx           = Wind.Velocity( 0, 1, y(iNd, iB), z(iNd, iB));
        phi(iNd, iB) = atan2(Vx, Vy(iNd));
    end
end

