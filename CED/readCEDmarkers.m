function [CEDStruct]  = readCEDmarkers(datapath)
% 
% readCEDevents.m reads all event data from one Spike2 file : events, level, 
% markers, textmarks, realmarks, wavemarks. It reads the time of each event,
% and the extended information, ie marker code, text, data, depending on 
% the event channel type.
% 
% Use as :
%       [CEDStruct]  = readCEDevents(datapath)
% 
% ## INPUT 
% datapath : path of the Spike2 file. It can be given with or without the 
%            extension .smr/.smrx. 
%
% ## OUTPUT
% CEDStruct : structure with all event timings from each
%             event/marker/realmark/textmark channel of the Spike2 file.
%
% Notes : 
% - The data file must not be opened in Spike2 when this script is running. 
% - readCEDevents.m needs the CEDS64ML interface library, and the Spike2 software.
%   The CEDS65ML library is automatically loaded in this script if the folder
%   'CED_library' is in the same path as this script. 
% - readCEDevents.m can only be ran on Windows (because of the CED library).
% - readCEDevents.m uses renamechan_CED.m. This function must stay in the same
%   folder as readCEDevents.m.
%
% paul.baudin@live.fr
%

% load CED library
if ~libisloaded('ceds64int')
    mfile_name  = fileparts(mfilename('fullpath'));
    cedpath     = fullfile(mfile_name, 'CED_library', 'CEDS64ML');
    if isfolder(cedpath)
        addpath(cedpath);
        CEDS64LoadLib(cedpath);
    else
        error('The CEDS64ML library must be loaded before calling this script');
    end
end

fprintf('Extracting events timings from  %s \n',datapath);

%remove extension if any
[folder, file]  = fileparts(datapath);
datapath        = fullfile(folder,file);

%find extension
temp     = dir([datapath, '.smr*']);%.smr* because some data are .smr and other .smrx
datapath = fullfile(temp.folder, temp.name);

%Open file
fid = CEDS64Open(datapath);
if fid<0
    error('Error while opening file : %s \nPossible issues : \n - the data path is wrong\n - the file is opened in Spike2.',datapath);
end

%recover "real" time
[~, timeinteger]    = CEDS64TimeDate(fid);
timems              = timeinteger(1)*10;
timesec             = timeinteger(2);
timemin             = timeinteger(3);
timehour            = timeinteger(4);
timeday             = timeinteger(5);
timemonth           = timeinteger(6);
timeyear            = timeinteger(7);
timestring          = sprintf('%d/%d/%d %d:%d:%d.%d',timeyear, timemonth, timeday, timehour, timemin, timesec, timems);
starttime           = datetime(timestring, 'Format', 'yyyy/MM/dd HH:mm:ss.SSS');
maxtime             = CEDS64TicksToSecs(fid, CEDS64MaxTime(fid));

CEDStruct.starttime = starttime;
CEDStruct.filename  = datapath;
CEDStruct.endtime   = starttime + seconds(maxtime);

%get channel nr
channr = CEDS64MaxChan(fid);

%go through each channel
i = 0;
for ichan = 1:channr
        
    [iType] = CEDS64ChanType(fid, ichan);
    if ~ismember(iType, 2:8) 
        continue
    end
    
    i = i + 1;
    
    %Get channel name and rename it if needed
    [~, name{i}] = CEDS64ChanTitle(fid, ichan);
    name{i}      = renamechan_CED(name{i}, ichan, name(1:i-1), true); %name{1:i-1} is [] if i = 1
    CEDStruct.markers.(name{i}).comment        = sprintf('Spike2 chan nr %d', ichan);
    
    %read timings
    [iread, events] = CEDS64ReadEvents(fid,ichan,100000,0);
    if iread > 100000
        error('Too many events in channel %d of file %s. \nIncrease the corresponding argument in CEDS64ReadEvents just above this line', ichan, datapath);
    end
    if iread < 0
        error('Error while loading events of channel %d in file %s', ichan, datapath);
    end
    fprintf('Found %d occurences of %s\n', iread, name{i});
    
    %store times
    if ~isempty(events)
        CEDStruct.markers.(name{i}).synctime = (CEDS64TicksToSecs(fid,events))';
        CEDStruct.markers.(name{i}).clock    = seconds(CEDStruct.markers.(name{i}).synctime) + CEDStruct.starttime;
    else
        CEDStruct.markers.(name{i}).synctime = [];
        CEDStruct.markers.(name{i}).clock = [];
    end
    
    % read extended infos depending on the channel type
    switch iType
        
        case 2
            CEDStruct.markers.(name{i}).chantype = 'event-';
            
        case 3
            CEDStruct.markers.(name{i}).chantype = 'event+';
            
        case 4 %level
            [~, marksinfos] = CEDS64ReadMarkers(fid, ichan,100000,0);
            if iread > 100000
                error('Too many events in channel %d of file %s. \nIncrease the corresponding argument in CEDS64ReadMarkers just above this line', ichan, datapath);
            end
            for ievent = 1:size(marksinfos, 1)
                CEDStruct.markers.(name{i}).level(ievent)    = marksinfos(ievent).m_Code1;
            end
            CEDStruct.markers.(name{i}).chantype = 'level';
            
        case 5 %markers
            [iread, marksinfos] = CEDS64ReadMarkers(fid, ichan,100000,0);
            if iread > 100000
                error('Too many events in channel %d of file %s. \nIncrease the corresponding argument in CEDS64ReadMarkers just above this line', ichan, datapath);
            end
            for ievent = 1:size(marksinfos, 1)
                CEDStruct.markers.(name{i}).code_1(ievent)   = marksinfos(ievent).m_Code1;
                CEDStruct.markers.(name{i}).code_2(ievent)   = marksinfos(ievent).m_Code2;
                CEDStruct.markers.(name{i}).code_3(ievent)   = marksinfos(ievent).m_Code3;
                CEDStruct.markers.(name{i}).code_4(ievent)   = marksinfos(ievent).m_Code4;
            end
            CEDStruct.markers.(name{i}).chantype = 'marker';
            
        case 6 %wavemarks
            [iread, marksinfos] = CEDS64ReadExtMarks(fid, ichan,100000,0);
            if iread > 100000
                error('Too many events in channel %d of file %s. \nIncrease the corresponding argument in CEDS64ReadExtMarks just above this line', ichan, datapath);
            end
            for ievent = 1:size(marksinfos, 1)
                CEDStruct.markers.(name{i}).code_1(ievent)   = marksinfos(ievent).m_Code1;
                CEDStruct.markers.(name{i}).code_2(ievent)   = marksinfos(ievent).m_Code2;
                CEDStruct.markers.(name{i}).code_3(ievent)   = marksinfos(ievent).m_Code3;
                CEDStruct.markers.(name{i}).code_4(ievent)   = marksinfos(ievent).m_Code4;
                CEDStruct.markers.(name{i}).data{ievent}     = marksinfos(ievent).m_Data;
            end
            CEDStruct.markers.(name{i}).chantype = 'wavemark';
            
        case 7 %realmarks
            [iread, marksinfos] = CEDS64ReadExtMarks(fid, ichan,100000,0);
            if iread > 100000
                error('Too many events in channel %d of file %s. \nIncrease the corresponding argument in CEDS64ReadExtMarks just above this line', ichan, datapath);
            end
            for ievent = 1:size(marksinfos, 1)
                CEDStruct.markers.(name{i}).code_1(ievent)   = marksinfos(ievent).m_Code1;
                CEDStruct.markers.(name{i}).code_2(ievent)   = marksinfos(ievent).m_Code2;
                CEDStruct.markers.(name{i}).code_3(ievent)   = marksinfos(ievent).m_Code3;
                CEDStruct.markers.(name{i}).code_4(ievent)   = marksinfos(ievent).m_Code4;
                CEDStruct.markers.(name{i}).data(ievent)     = marksinfos(ievent).m_Data;
            end
            CEDStruct.markers.(name{i}).chantype = 'realmark';
            
        case 8 %textmarks
            [iread, marksinfos] = CEDS64ReadExtMarks(fid, ichan,100000,0);
            if iread > 100000
                error('Too many events in channel %d of file %s. \nIncrease the corresponding argument in CEDS64ReadExtMarks just above this line', ichan, datapath);
            end
            for ievent = 1:size(marksinfos, 1)
                CEDStruct.markers.(name{i}).code_1(ievent)   = marksinfos(ievent).m_Code1;
                CEDStruct.markers.(name{i}).code_2(ievent)   = marksinfos(ievent).m_Code2;
                CEDStruct.markers.(name{i}).code_3(ievent)   = marksinfos(ievent).m_Code3;
                CEDStruct.markers.(name{i}).code_4(ievent)   = marksinfos(ievent).m_Code4;
                CEDStruct.markers.(name{i}).text{ievent}     = marksinfos(ievent).m_Data;
            end
            CEDStruct.markers.(name{i}).chantype = 'textmark';
    end
        
end %ichan
    
if i == 0
    error('no event or marker channels found in %s', datapath);
else
    fprintf('Found %d event/marker channels in %s\n', i, datapath);
end

CEDStruct.markers = orderfields(CEDStruct.markers);

CEDS64Close(fid);