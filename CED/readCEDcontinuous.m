function chandata = readCEDcontinuous(datapath, channame)
% 
% readCEDwaveforms.m reads the Spike2 channel indicated in input, if it is a
% waveform channel or a realwave channel. It outputs the data in a Fieldtrip's
% raw format.
% 
% Use as :
%       chandata = readCEDwaveforms(datapath, channame)
% 
% ## INPUT 
% datapath : path of the Spike2 file. It can be given with or without the 
%            extension .smr/.smrx. 
% channame : name of the Spike2 channel to load. This name is sensible to 
%            upper or lower cases. It can contain white space or special  
%            characters. If several channels have the same name, only 
%            the last one (with the higher Spike2 chan nr) is loaded
%
% ## OUTPUT
% chandata : data in a Fieldtrip's raw format.
%
% Notes : 
% - The data file must not be opened in Spike2 when this script is running. 
% - readCEDwaveforms.m needs the CEDS64ML interface library, and the Spike2 software.
%   The CEDS65ML library is automatically loaded in this script if the folder
%   'CED_library' is in the same path as this script.
% - readCEDwaveforms.m can only be ran on Windows (because of the CED library).
% 
% paul.baudin@live.fr
% 

% load CED library
if ~libisloaded('ceds64int')
    mfile_name = fileparts(mfilename('fullpath'));
    cedpath    = fullfile(mfile_name, 'CED_library', 'CEDS64ML');
    if isfolder(cedpath)
        addpath(cedpath);
        CEDS64LoadLib(cedpath);
    else
        error('The CEDS64ML library must be loaded before calling this script');
    end
end

chandata        = [];
n_loaded_chan   = 0;

%remove extension if any
[folder, name]  = fileparts(datapath);
datapath        = fullfile(folder,name);

%find extension
temp = dir([datapath, '.smr*']);%.smr* because some data are .smr and other .smrx
if isempty(temp)
    error('Cannot find data in : %s', datapath);
end
datapath = fullfile(temp.folder, temp.name);
fprintf('Loading data from %s \nchannel %s\n', datapath, channame);

%open the file
fid = CEDS64Open(datapath);
if fid<0
    error('Error while opening file : %s \nPossible issues : \n - the data path is wrong\n - the file is opened in Spike2.',datapath);
end

% Go through each Spike2 channel
for ichan = 1:CEDS64MaxChan(fid)
    
    [~, chantitle] = CEDS64ChanTitle(fid, ichan);
    
    if strcmp(chantitle, channame)
        
        [iType] = CEDS64ChanType(fid, ichan);
        
        if ismember(iType, [1, 9]) %Waveform or RealWave
            chandata = [];
            n_loaded_chan = n_loaded_chan+1;
            
            %load data
            chanFs              = 1 / (CEDS64ChanDiv(fid, ichan) * CEDS64TicksToSecs(fid,1)); %inverse de : nb de tick par sample * durée du tick
            maxtime             = CEDS64TicksToSecs(fid, CEDS64MaxTime(fid));
            n_samples           = round(maxtime*chanFs);
            
            [n_samples_read, chandata.trial{1}, first_time] = CEDS64ReadWaveF(fid, ichan, n_samples, 0);
            first_time = CEDS64TicksToSecs(fid, first_time);
            
            %fill data if it does not start at zero, to make it easier to
            %deal with samples during the analysis
            if first_time > 0 + 1/chanFs
                toadd             = nan(size(first_time:-1/chanFs:0));
                n_samples_read    = n_samples_read + length(toadd);
                chandata.trial{1} = [toadd'; chandata.trial{1}];
                first_time        = 0;
            end
            
            %reconstruct discontinuous data, because CEDS64ReadWaveF stops
            %readind data when there is a discontinuity
            while n_samples_read < n_samples
                new_start = (n_samples_read + 1) / chanFs + first_time;
                [n_samples_read2, trialtemp, first_time2] = CEDS64ReadWaveF(fid, ichan, n_samples, CEDS64SecsToTicks(fid, new_start));
                
                if n_samples_read2 == 0
                    break
                end
                
                first_time2         = CEDS64TicksToSecs(fid, first_time2);
                missingtimes        = new_start : 1/chanFs : first_time2 - 1/chanFs;
                chandata.trial{1}   = [chandata.trial{1}; nan(size(missingtimes))';  trialtemp];
                n_samples_read      = n_samples_read + n_samples_read2 + length(missingtimes);
            end
                
            chandata.trial{1}        = chandata.trial{1}';
            chandata.fsample         = chanFs;
            chandata.sampleinfo      = [1,n_samples];
            
            chandata.time{1}         = first_time + (0 : 1/chanFs : (length(chandata.trial{1})-1)/chanFs); %-1 because begins at zero
            chandata.label{1}        = chantitle;
            [~,chandata.chanunit{1}] = CEDS64ChanUnits(fid,ichan);
            
        else
            
            if iType == 0
                continue
            end
            error('Channel %s is not a waveform or a realwave', channame);
            
        end
                
    end %strcmp channame chantitle
    
end %ichan

if n_loaded_chan == 0
    warning('Channel %s was not found in %s', channame, datapath);
elseif n_loaded_chan > 1
    warning('%d channels have the name %s. Only the last one (with the higher Spike2 chan nr) was loaded', n_loaded_chan, channame);
end

CEDS64Close(fid);