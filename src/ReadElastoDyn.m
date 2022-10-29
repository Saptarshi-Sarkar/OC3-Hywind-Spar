function [ElastoDyn] = ReadElastoDyn()

% Read elastodyn file provided by NREL for all the inputs

fid = fopen('./5MW_Baseline/NRELOffshrBsline5MW_OC3Hywind_ElastoDyn.dat','r');
Info = textscan(fid,'%f %s',17,'delimiter','\n','HeaderLines',27);
% Initial conditions
ElastoDyn.OoPDefl = Info{1,1}(1);
ElastoDyn.IPDefl  = Info{1,1}(2);
ElastoDyn.BlPitch(1) = Info{1,1}(3)/180*pi;
ElastoDyn.BlPitch(2) = Info{1,1}(4)/180*pi;
ElastoDyn.BlPitch(3) = Info{1,1}(5)/180*pi;
ElastoDyn.Azimuth = Info{1,1}(7);
ElastoDyn.RotSpeed = Info{1,1}(8)/60*2*pi;                % in rad/sec
ElastoDyn.NacYaw   = Info{1,1}(9);
ElastoDyn.TTDspFA  = Info{1,1}(10);
ElastoDyn.TTDspSS  = Info{1,1}(11);
ElastoDyn.PtfmSurge = Info{1,1}(12);
ElastoDyn.PtfmSway = Info{1,1}(13);
ElastoDyn.PtfmHeave = Info{1,1}(14);
ElastoDyn.PtfmRoll = Info{1,1}(15);
ElastoDyn.PtfmPitch = Info{1,1}(16);
ElastoDyn.PtfmYaw = Info{1,1}(17);

Info = textscan(fid,'%f %s',26,'delimiter','\n','HeaderLines',1);
% Turbine Configuration 
ElastoDyn.NumBl = Info{1,1}(1);
ElastoDyn.TipRad = Info{1,1}(2);
ElastoDyn.HubRad = Info{1,1}(3);
ElastoDyn.PreCone(1) = Info{1,1}(4)/180*pi;
ElastoDyn.PreCone(2) = Info{1,1}(5)/180*pi;
ElastoDyn.PreCone(3) = Info{1,1}(6)/180*pi;
ElastoDyn.HubCM = Info{1,1}(7);
ElastoDyn.OverHang = Info{1,1}(11);
ElastoDyn.ShftGagL = Info{1,1}(12);
ElastoDyn.ShftTilt = Info{1,1}(13)/180*pi;
ElastoDyn.NacCMxn = Info{1,1}(14);
ElastoDyn.NacCMyn = Info{1,1}(15);
ElastoDyn.NacCMzn = Info{1,1}(16);
ElastoDyn.NcIMUxn = Info{1,1}(17);
ElastoDyn.NcIMUyn = Info{1,1}(18);
ElastoDyn.NcIMUzn = Info{1,1}(19);
ElastoDyn.Twr2Shft = Info{1,1}(20);
ElastoDyn.TwrHt = Info{1,1}(21);
ElastoDyn.TowerBsHt = Info{1,1}(22);  
ElastoDyn.PtfmCMxt = Info{1,1}(23);    
ElastoDyn.PtfmCMyt = Info{1,1}(24);   
ElastoDyn.PtfmCMzt = Info{1,1}(25);   
ElastoDyn.PtfmRefzt = Info{1,1}(26);   
Info = textscan(fid,'%f %s',13,'delimiter','\n','HeaderLines',1);
% Mass and Inertia
ElastoDyn.HubMass = Info{1,1}(4);
ElastoDyn.HubIner = Info{1,1}(5);
ElastoDyn.GenIner = Info{1,1}(6);
ElastoDyn.NacMass = Info{1,1}(7);
ElastoDyn.NacYIner = Info{1,1}(8);
ElastoDyn.YawBrMass = Info{1,1}(9);
ElastoDyn.PtfmMass = Info{1,1}(10);
ElastoDyn.PtfmRIner = Info{1,1}(11);
ElastoDyn.PtfmPIner = Info{1,1}(12);
ElastoDyn.PtfmYIner = Info{1,1}(13);

% Drivetrain
Info = textscan(fid,'%f %s',4,'delimiter','\n','HeaderLines',15);
ElastoDyn.GBoxEff = Info{1,1}(1);
ElastoDyn.GBRatio = Info{1,1}(2);
ElastoDyn.DTTorSpr = Info{1,1}(3);
ElastoDyn.DTTorDmp = Info{1,1}(4);

fclose(fid);





