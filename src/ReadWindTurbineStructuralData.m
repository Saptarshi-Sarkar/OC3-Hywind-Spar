function [Blade, Tower] = ReadWindTurbineStructuralData()

% Read structural data files provided by NREL

% Read blade data 
fid = fopen('./5MW_Baseline/NRELOffshrBsline5MW_Blade.dat','r');
Info = textscan(fid,'%f %s',4,'delimiter','\n','HeaderLines',3);
Blade.NoOfPoints = Info{1,1}(1);
Blade.DampBF1 = Info{1,1}(2);
Blade.DampBF2 = Info{1,1}(3);
Blade.DampBE1 = Info{1,1}(4);
formatSpec = repmat('%f', 1, 6);
Data = textscan(fid,formatSpec,Blade.NoOfPoints,'delimiter','\n','HeaderLines',9);
Blade.Data = cell2mat(Data);
Blade.Data(:,3) = Blade.Data(:,3)*pi/180;
ModeData = textscan(fid,'%f %s',15,'delimiter','\n','HeaderLines',1);
Blade.FlapMode1 = [flipud(ModeData{1,1}(1:5));0;0]';
Blade.FlapMode2 = [flipud(ModeData{1,1}(6:10));0;0]'; 
Blade.EdgeMode1 = [flipud(ModeData{1,1}(11:15));0;0]';
fclose(fid);

% Read Tower data
fid = fopen('./5MW_Baseline/NRELOffshrBsline5MW_OC3Hywind_ElastoDyn_Tower.dat','r');
Info = textscan(fid,'%f %s',5,'delimiter','\n','HeaderLines',3);
Tower.NoOfPoints = Info{1,1}(1);
Tower.Damp_TFA1 = Info{1,1}(2);
Tower.Damp_TFA2 = Info{1,1}(3);
Tower.Damp_TSS1 = Info{1,1}(4);
Tower.Damp_TSS2 = Info{1,1}(5);
formatSpec = repmat('%f', 1, 4);
Data = textscan(fid,formatSpec,Tower.NoOfPoints,'delimiter','\n','HeaderLines',11);
Tower.Data = cell2mat(Data);
ModeData = textscan(fid,'%f %s',10,'delimiter','\n','HeaderLines',1);
Tower.FAMode1 = [flipud(ModeData{1,1}(1:5));0;0]';
Tower.FAMode2 = [flipud(ModeData{1,1}(6:10));0;0]';
ModeData = textscan(fid,'%f %s',10,'delimiter','\n','HeaderLines',1);
Tower.SSMode1 = [flipud(ModeData{1,1}(1:5));0;0]';
Tower.SSMode2 = [flipud(ModeData{1,1}(6:10));0;0]';
fclose(fid);



