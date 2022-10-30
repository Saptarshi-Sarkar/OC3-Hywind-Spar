function CheckInterpPoints(GridZ, GridY, ZBlNode_Y, ZBlNode_Z)
% Wind.z = vertical
% Wind.y = horizontal
% ZBlNode_Y = vertical
% ZBlNode_Z = horizontal

for iBd = 1:3
    if max(ZBlNode_Y(:,iBd)) > GridY(end) 
        error('Highest vertical interpolation point %f beyond grid. Increase grid height', max(ZBlNode_Y(:,iBd)));
    elseif min(ZBlNode_Y(:,iBd)) < GridY(1)
        error('Lowest vertical interpolation point %f beyond grid. Increase grid height', min(ZBlNode_Y(:,iBd)));        
    elseif max(ZBlNode_Z(:,iBd)) > GridZ(end) 
        error('Horizontal interpolation point beyond %f grid. Increase grid width', max(ZBlNode_Z(:,iBd)));
    elseif min(ZBlNode_Z(:,iBd)) < GridZ(1)
        error('Horizontal interpolation point beyond %f grid. Increase grid width', min(ZBlNode_Z(:,iBd)));
    end
end
end

