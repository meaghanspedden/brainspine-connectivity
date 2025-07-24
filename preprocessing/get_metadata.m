% metadata
function [runs, ref,heartidx,badsenslabels] = get_metadata(subID,analysis)

runs=struct();

badsenslabels = {};


if strcmp(subID, 'OP00212')
    if strcmp(analysis,'static')
        runs.opm={'001', '002', '003', '004', '005'};
        runs.emg={'01', '04', '06', '08', '10'};
    else
        runs.opm={'002', '003', '004', '005', '006'};
        runs.emg={'03', '05', '07', '09', '11'};
    end
    ref='25';
    heartidx=2;


elseif strcmp(subID, 'OP00213')
    if strcmp(analysis, 'static')
        runs.opm={'001', '002', '003', '004', '005'};
        runs.emg={'01', '03', '05', '07', '09'};
    else
        runs.opm={'001', '002', '003', '004', '005'};
        runs.emg={'02', '04', '06', '08', '10'};
    end
    ref='25';
    heartidx=2;
    badsenslabels = {'28'};  % identified post hoc


elseif strcmp(subID, 'OP00214')

    if strcmp(analysis, 'static')
        runs.opm={'001', '002', '003', '004', '005'};
        runs.emg={'01', '03', '05', '07', '09'};
    else
        runs.opm={'001', '002', '003', '004', '005'};
        runs.emg={'02', '04', '06', '08', '10'};
    end
    ref=[];%'27'; it was 27 but dropped out...
    heartidx=2;
    badsenslabels = {'27'};


elseif strcmp(subID, 'OP00215')
    if strcmp(analysis, 'static')
        runs.opm={'001', '002', '003', '004'};
        runs.emg={'01', '03', '05', '07'};
    else
        runs.opm={'001', '002', '003', '004'};
        runs.emg={'02', '04', '06', '08'};
    end
    ref='37';
    heartidx=2;
    badsenslabels = {'42'};


elseif strcmp(subID, 'OP00219')

    if strcmp(analysis, 'static')
        runs.opm={'001', '002', '003', '004'};
        runs.emg={'1', '3', '5', '7'};
    else
        runs.opm={'001', '002', '003', '004'};
        runs.emg={'2', '4', '6', '8'};
    end

    ref='26';
    heartidx=3;


elseif strcmp(subID, 'OP00220')

    if strcmp(analysis, 'static')
        runs.opm={'001', '002', '003', '004'};
        runs.emg={'1', '3', '5', '7'};
    else
        runs.opm={'001', '002', '003', '004'};
        runs.emg={'2', '4', '6', '8'};
    end
    ref='26';
    heartidx=3;

    badsenslabels = {'X27'};


elseif strcmp(subID, 'OP00221')
    if strcmp(analysis, 'static')
        runs.opm={'001', '002', '003', '004'};
        runs.emg={'1', '3', '5', '7'};
    else
        runs.opm={'001', '002', '003', '004'};
        runs.emg={'2', '4', '6', '8'};
    end
    ref='26';
    heartidx=3;
    badsenslabels = {'X20', 'X46'};


end