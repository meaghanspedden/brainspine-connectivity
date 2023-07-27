function [trigtimes, EMGstruct]= readSpikeData(datapath,channame)

if isempty(getenv('CEDS64ML'))
    setenv('CEDS64ML', 'C:\path\to\your\folder\CEDS64ML');
end
cedpath = getenv('CEDS64ML');
addpath(cedpath);
% load ceds64int.dll


CEDS64LoadLib( cedpath );
EMGstruct = readCEDcontinuous(datapath, channame);


[CEDStruct]  = readCEDmarkers(datapath);
try
trigtimes=CEDStruct.markers.Matlab_tr.synctime;
catch
trigtimes=CEDStruct.markers.trig.synctime;

end