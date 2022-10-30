function [Cl, Cd] = LiftDragCoeffInterp(AOA, FoilNo, Airfoils)
%#Codegen
if FoilNo == 1
    Cl = interp1(Airfoils.Cylinder1(:,1), Airfoils.Cylinder1(:,2), AOA,'linear','extrap');
    Cd = interp1(Airfoils.Cylinder1(:,1), Airfoils.Cylinder1(:,3), AOA,'linear','extrap');
elseif FoilNo == 2
    Cl = interp1(Airfoils.Cylinder2(:,1), Airfoils.Cylinder2(:,2), AOA,'linear','extrap');
    Cd = interp1(Airfoils.Cylinder2(:,1), Airfoils.Cylinder2(:,3), AOA,'linear','extrap');
elseif FoilNo == 3
    Cl = interp1(Airfoils.DU21_A17(:,1), Airfoils.DU21_A17(:,2), AOA,'linear','extrap');
    Cd = interp1(Airfoils.DU21_A17(:,1), Airfoils.DU21_A17(:,3), AOA,'linear','extrap');
elseif FoilNo == 4
    Cl = interp1(Airfoils.DU25_A17(:,1), Airfoils.DU25_A17(:,2), AOA,'linear','extrap');
    Cd = interp1(Airfoils.DU25_A17(:,1), Airfoils.DU25_A17(:,3), AOA,'linear','extrap');
elseif FoilNo == 5
    Cl = interp1(Airfoils.DU30_A17(:,1), Airfoils.DU30_A17(:,2), AOA,'linear','extrap');
    Cd = interp1(Airfoils.DU30_A17(:,1), Airfoils.DU30_A17(:,3), AOA,'linear','extrap');
elseif FoilNo == 6
    Cl = interp1(Airfoils.DU35_A17(:,1), Airfoils.DU35_A17(:,2), AOA,'linear','extrap');
    Cd = interp1(Airfoils.DU35_A17(:,1), Airfoils.DU35_A17(:,3), AOA,'linear','extrap');
elseif FoilNo == 7
    Cl = interp1(Airfoils.DU40_A17(:,1), Airfoils.DU40_A17(:,2), AOA,'linear','extrap');
    Cd = interp1(Airfoils.DU40_A17(:,1), Airfoils.DU40_A17(:,3), AOA,'linear','extrap');
elseif FoilNo == 8
    Cl = interp1(Airfoils.NACA64_A17(:,1), Airfoils.NACA64_A17(:,2), AOA,'linear','extrap');
    Cd = interp1(Airfoils.NACA64_A17(:,1), Airfoils.NACA64_A17(:,3), AOA,'linear','extrap');
else
    Cl = 0;
    Cd = 0;
end

end

