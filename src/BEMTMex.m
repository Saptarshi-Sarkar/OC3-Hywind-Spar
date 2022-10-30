function [ px, py, Mz, phi] = BEMTMex(q_GeAz, VMB, UMB, phi, nb_Load, BlSpn, FlexBlSpn, TipRad, HubRad, Cord, Twist,...
                                    Pitch, chi0, FoilNo, Airfoils, AeroCentJ1, AeroCentJ2)
    rho = 1.225;
    
    px = zeros(nb_Load,3);
    py = zeros(nb_Load,3);
    Mz = zeros(nb_Load,3);
    
    for B = 1:3
        Azimuth = q_GeAz + (B-1)*2*pi/3;
        
        for iNd = 2:1:nb_Load-1
            
            Solid = (3/2/pi)*Cord(iNd)/FlexBlSpn(iNd, B);
            Vx = VMB(iNd,1,B) - UMB(iNd,1,B);
            Vy = VMB(iNd,2,B) + UMB(iNd,2,B);
            
            theta = Pitch(B) + Twist(iNd);
                           
            Residue =  CallStateResidual(phi(iNd, B), theta, Vx, Vy, FoilNo(iNd), Airfoils, TipRad, HubRad, BlSpn(iNd), Solid, Azimuth, chi0);                
            if Vx == 0 || Vy == 0
                phi(iNd, B) = atan2(Vx, Vy);
                [a, aa, Cx, Cy, Cl, Cd] = CalcOutput(phi(iNd, B), theta, FoilNo(iNd), Airfoils, TipRad, HubRad, BlSpn(iNd), Solid, Azimuth, chi0);
            elseif abs(Residue) < 1e-3
                [a, aa, Cx, Cy, Cl, Cd] = CalcOutput(phi(iNd, B), theta, FoilNo(iNd), Airfoils, TipRad, HubRad, BlSpn(iNd), Solid, Azimuth, chi0);
            else
                phi(iNd, B) =  UpdateInflowAngle(theta, Vx, Vy, FoilNo(iNd), Airfoils, TipRad, HubRad, BlSpn(iNd), Solid, Azimuth, chi0);
                [a, aa, Cx, Cy, Cl, Cd] = CalcOutput(phi(iNd, B), theta, FoilNo(iNd), Airfoils, TipRad, HubRad, BlSpn(iNd), Solid, Azimuth, chi0);
            end           
            
            % Final Velocities and forces %            
            W = (Vx*(1-a))^2 + (Vy*(1+aa))^2;
            Const  = 0.5*rho*Cord(iNd)*W;

            px(iNd,B) = Cx*Const;
            py(iNd,B) = -Cy*Const;
            
            Mz(iNd,B) = Cl*Const*AeroCentJ2(iNd) - Cd*Const*AeroCentJ1(iNd); % + Cm(iNd)*Const*Cord(iNd)           
        end
            
    end
end

% Force coefficients %
function [Cx, Cy] = ForceCoefficents(AOA, phi, FoilNo, Airfoils)
    [Cl, Cd] = LiftDragCoeffInterp(AOA, FoilNo, Airfoils);
    Cx = Cl*cos(phi) + Cd*sin(phi);
    Cy = Cl*sin(phi) - Cd*cos(phi);
end

% Hub/Tip Losses %
function F = HubTipLoss(phi, TipRad, HubRad, BlSpn)
    abssinphi = abs(sin(phi));
    ftip = 1.5*(TipRad - BlSpn)/(BlSpn*abssinphi);
    Ftip = 2.0/pi*acos(min(1,exp(-ftip)));

    fhub = 1.5*(BlSpn - HubRad)/(HubRad*abssinphi);
    Fhub = 2.0/pi*acos(min(1,exp(-fhub)));
    F    = Ftip*Fhub;
end

% Residue %
function Residue =  CallStateResidual(phi, theta, Vx, Vy, FoilNo, Airfoils, TipRad, HubRad, BlSpn, Solid, Azimuth, chi0)
    AOA      = phi - theta;
    [Cx, Cy] = ForceCoefficents(AOA, phi, FoilNo, Airfoils);
    F        = HubTipLoss(phi, TipRad, HubRad, BlSpn);

    cosphi = cos(phi);
    sinphi = sin(phi);

    % Non dimensional parameters %
    if sinphi == 0
        k = realmax;
    else
        k = Solid*Cx/(4*F*sinphi^2);
    end
    a = 0;
%     aa = 1;
    % Different equation depending on solution region %
    if phi > 0
        if k <= 2/3 % momentum region
            if k == -1
                k = k - 0.1;
            end
            a = k/(1.0+k);
            a = max(a, -10.0);          % Patch
        else     % emperical region
            g1 = 2.0*F*k - (10/9 - F);
            g2 = 2.0*F*k - F*(4/3 - F);
            g3 = 2.0*F*k - (25/9 - 2.0*F);
            if (abs(g3) < 1e-6)
                a = 1 - 1/(2*sqrt(g2));
            else
                a = (g1 - sqrt(g2))/g3;
            end
        end
    elseif phi < 0 
        if k > 1          % propeller brake region
        a = k/(k-1.0);
        a = min(a,10);             % Patch
        else
        a = 0.0;
        end
    elseif phi == 0
        a  = 0;
        aa = 0;
        Residue = sinphi/(1-a) - Vx/Vy*cosphi/(1+aa);
        return
    end

    %Pitt and Peters yaw correction model %
    x = (0.6*a+1)*(chi0);
    a = a*(1.0 + 15*pi/64*tan(x/2.0)*BlSpn/TipRad*sin(Azimuth));

    % Tangential induction factor %
    if sinphi*cosphi == 0
        kk = realmax;
    else
        kk = Solid*Cy/(4*F*sinphi*cosphi);
    end
    
    if kk == 1.0
        a = 0.0;
        aa = 0.0;
    else
        aa = kk/(1-kk);
        if abs(aa) > 10.0
            aa = 10*sign(aa); % Patch
        end
    end
    
    %% Residue %
    if a == 1
        Residue = - Vx/Vy*cosphi/(1+aa);
    else
        Residue = sinphi/(1-a) - Vx/Vy*cosphi/(1+aa);
    end
end

% Update Inflow Angle %
function phi = UpdateInflowAngle(theta, Vx, Vy, FoilNo, Airfoils, TipRad, HubRad, BlSpn, Solid, Azimuth, chi0)
    options = optimset('Display','notify','TolFun',1e-3);
    ep = 1e-6;
    Res = @(phi) CallStateResidual(phi, theta, Vx, Vy, FoilNo, Airfoils, TipRad, HubRad, BlSpn, Solid, Azimuth, chi0);

    if Res(pi/2)*Res(ep) < 0
        x0 = ep;      
        x1 = pi/2;    
        phi =  fzero(Res,[x0, x1], options);  % Standard solution region
    elseif Res(-pi/4)*Res(-ep) < 0
        x0 = -pi/4;
        x1 = -ep;
        phi =  fzero(Res,[x0, x1], options);  % Propeller brake region
    else
        disp(Vx);
        disp(Vy);
        x0 = pi/2;
        x1 = pi-ep;
        phi =  fzero(Res,[x0, x1], options);  % Last solution region
    end

end

% Calculate Output %
function [a, aa, Cx, Cy, Cl, Cd] =  CalcOutput(phi, theta, FoilNo, Airfoils, TipRad, HubRad, BlSpn, Solid, Azimuth, chi0)
    AOA      = phi - theta;
    [Cl, Cd] = LiftDragCoeffInterp(AOA, FoilNo, Airfoils);
    [Cx, Cy] = ForceCoefficents(AOA, phi, FoilNo, Airfoils);
    F        = HubTipLoss(phi, TipRad, HubRad, BlSpn);

    cosphi = cos(phi);
    sinphi = sin(phi);

    % Non dimensional parameters %
    if sinphi == 0
        k = realmax;
    else
        k = Solid*Cx/(4*F*sinphi^2);
    end
    
    % Different equation depending on solution region %
    if phi > 0
        if k <= 2/3 % momentum region
            if k == -1
                k = k - 0.1;
            end
            a = k/(1.0+k);
            a = max(a, -10.0);          % Patch
        else     % emperical region
            g1 = 2.0*F*k - (10/9 - F);
            g2 = 2.0*F*k - F*(4/3 - F);
            g3 = 2.0*F*k - (25/9 - 2.0*F);
            if (abs(g3) < 1e-6)
                a = 1 - 1/(2*sqrt(g2));
            else
                a = (g1 - sqrt(g2))/g3;
            end
        end
    elseif phi < 0 
        if k > 1          % propeller brake region
        a = k/(k-1.0);
        a = min(a,10);             % Patch
        else
        a = 0.0;
        end
    else
        a  = 0;
        aa = 0;
        return
    end

    %Pitt and Peters yaw correction model %
    x = (0.6*a+1)*(chi0);
    a = a*(1.0 + 15*pi/64*tan(x/2.0)*BlSpn/TipRad*sin(Azimuth));

    % Tangential induction factor %
    if sinphi*cosphi == 0
        kk = realmax;
    else
        kk = Solid*Cy/(4*F*sinphi*cosphi);
    end

    if kk == 1.0
        a = 0.0;
        aa = 0.0;
    else
        aa = kk/(1-kk);
        if abs(aa) > 10.0
            aa = 10*sign(aa); % Patch
        end
    end

end