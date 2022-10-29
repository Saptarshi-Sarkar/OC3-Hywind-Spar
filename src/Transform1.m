function [ Out ] = Transform1( RotMat, BasisVec )
%#Codegen
% RotMat : Rotation  matrix (3x3)
% BasisVec : Matrix with basis vectors, each row containing a basis vector
Out = zeros(3,3);
Out(1,:) = RotMat(1,1)*BasisVec(1,:) + RotMat(1,2)*BasisVec(2,:) + RotMat(1,3)*BasisVec(3,:);
Out(2,:) = RotMat(2,1)*BasisVec(1,:) + RotMat(2,2)*BasisVec(2,:) + RotMat(2,3)*BasisVec(3,:);
Out(3,:) = RotMat(3,1)*BasisVec(1,:) + RotMat(3,2)*BasisVec(2,:) + RotMat(3,3)*BasisVec(3,:);

end

