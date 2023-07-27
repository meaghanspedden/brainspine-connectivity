% This function opens a file and reads the data from a waveform channel as 32-bit floats. Then downsamples the data and saves it to a new file. It does not alter the orignal data.
% clear workspace

file2open='fingertapping_00.smrx';


clear;
% add path to CED code
if isempty(getenv('CEDS64ML'))
    setenv('CEDS64ML', 'C:\path\to\your\folder\CEDS64ML');
end
cedpath = getenv('CEDS64ML');
addpath(cedpath);
% load ceds64int.dll


CEDS64LoadLib( cedpath );
% Open a file
fhand1 = CEDS64Open(file2open);
if (fhand1 <= 0); unloadlibrary ceds64int; return; end


ch=1;
maxTimeTicks = CEDS64ChanMaxTime( fhand1, ch )+1; % +1 so the read gets the last point as reads up to

[ fRead, fVals4, fTime ] = CEDS64ReadWaveF( fhand1, ch, 200000000, 0, maxTimeTicks );


rate = CEDS64IdealRate( fhand1, 1 ); %samp freq


%need to read trigger and clip
[ fRead, fVals4] = CEDS64ReadMarkers( fhand1, ch, 200000000, 0, -1);





% close all the files
CEDS64CloseAll();
% unload ceds64int.dll
unloadlibrary ceds64int;