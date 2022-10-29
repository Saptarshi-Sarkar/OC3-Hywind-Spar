function [qd, Controls] = RHS(t, q, Controls, DOFs, ElastoDyn, Airfoils, Twr, Bld, Platform, Wind, WindNom, wave, mooring_load_ptr, Servo)

nDOFs = DOFs.nDOFs;
[C,f,Controls] = SystemMatrices(t, q, Controls, DOFs, ElastoDyn, Airfoils, Twr, Bld, Platform, Wind, WindNom, wave, mooring_load_ptr, Servo);
qd(1,1:nDOFs) = q(1,nDOFs+1:2*nDOFs);
qd(1,nDOFs+1:2*nDOFs) = C\f;

end

