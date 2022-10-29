function [ Out ] = Transform2( RotMat, BasisVec )
%#Codegen
% RotMat : Rotation  matrix (3 x 3 x N)
% BasisVec : Matrix (3 x 3) with basis vectors, each row containing a basis vector
Out = zeros(size(RotMat));
for i = 1:size(RotMat,3)
    Out(1,:,i) = RotMat(1,1,i)*BasisVec(1,:) + RotMat(1,2,i)*BasisVec(2,:) + RotMat(1,3,i)*BasisVec(3,:);
    Out(2,:,i) = RotMat(2,1,i)*BasisVec(1,:) + RotMat(2,2,i)*BasisVec(2,:) + RotMat(2,3,i)*BasisVec(3,:);
    Out(3,:,i) = RotMat(3,1,i)*BasisVec(1,:) + RotMat(3,2,i)*BasisVec(2,:) + RotMat(3,3,i)*BasisVec(3,:);
end
end