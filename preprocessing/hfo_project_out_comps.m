function [mfD,Yinds] = hfo_project_out_comps(S)
% Project out any components of interest using an offline correction.
% FORMAT D = hfo_project_out_comps(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                                - Default: no Default
%   S.X             - Matrix of components to project out (nchans x
%                                           projectors)              - Default: REQUIRED                                                                  
%   S.channels      - indicies of the channels for correction        - Default: empty
%   S.chunkSize     - max memory usage(for large datasets)           - Default 512(MB)
%   S.balance       - logical to update forward model                - Default 1
%   S.prefix        - prefix to filename                             - Default 'p'
% Output:
%   D               - denoised MEEG object (also written to disk)
%   Yinds           - the indices of filtered channels
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging

%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),             error(errorMsg); end
if ~isfield(S, 'X'),             error('please supply deisgn matrix'); end
if ~isfield(S, 'channels'),      S.channels = []; end
if ~isfield(S, 'chunkSize'),     S.chunkSize = 512; end
if ~isfield(S, 'balance'),       S.balance = 1; end
if ~isfield(S, 'prefix'),        S.prefix = 'p'; end

%-Prep sensor identification between SPM and FIELDTRIP
%--------------------------------------------------------------------------

s = sensors(S.D,'MEG');
if isempty(s)==1
    error('Could not find sensor positions')
end

if isempty(S.channels)
    S.channels = setxor(S.D.selectchannels('MEG'),S.D.badchannels);
end

spm_labs = S.D.chanlabels(S.channels);
[ind, ft_labs] = match_str(spm_labs,s.label);


%-Get design matrix
%--------------------------------------------------------------------------
X = S.X;

%-Compute projector
%--------------------------------------------------------------------------
M = eye(size(X,1))-X*pinv(X);

%-Get Data indices
%--------------------------------------------------------------------------
Yinds = S.channels;

if numel(Yinds)~=size(X,1)
    error('data size ~= number sensors in components');
end

%-create ouput dataset object
%--------------------------------------------------------------------------
fprintf('Creating output dataset\n');
outname = fullfile(path(S.D),[S.prefix fname(S.D)]);
mfD = clone(S.D,outname);
mfD.save();

%- Work out chunk size
%--------------------------------------------------------------------------
chunkSamples= round(S.chunkSize/(8*size(S.D,1))*1e6);
begs=1:chunkSamples:size(S.D,2);
ends = (begs+chunkSamples-1);
if(ends(end)>size(S.D,2))
    ends(end)= size(S.D,2);
end

%-Run on channels needing correction
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Processing Data',spm('time'));

for j=1:size(S.D,3)
    mk= S.D(Yinds,1,j);
    sk = zeros(length(Yinds),1);
    count=1;
    
    for i =1:length(begs)
        inds = begs(i):ends(i);
        
        % mfD(Yinds,inds,j)=M*S.D(Yinds,inds,j) is slow (disk read)
        out = S.D(:,inds,j);
        Y=out(Yinds,:);
        out(Yinds,:)=M*Y;
        mfD(:,inds,j)=out;
        
    end
    
end

%-Update forward modelling information
%--------------------------------------------------------------------------
if (S.balance)
    fprintf('%-40s: %30s\n','Updating Sensor Information',spm('time'));
    grad = mfD.sensors('MEG');
    tmpTra= eye(size(grad.chanpos,1));
    tmpTra(ft_labs,ft_labs)=M;
    grad.tra                = tmpTra*grad.tra;
    grad.balance.previous   = grad.balance.current;
    grad.balance.current    = 'reject_comps';
    mfD = sensors(mfD,'MEG',grad);
    % Check if any information in D.inv needs updating.
    % TODO: Update to support multiple invs/forwards/modalities
    if isfield(mfD,'inv')
        if isfield(mfD.inv{1},'gainmat')
            fprintf(['Clearing current forward model, please recalculate '...
                'with spm_eeg_lgainmat\n']);
            mfD.inv{1} = rmfield(mfD.inv{1},'gainmat');
        end
        if isfield(mfD.inv{1},'datareg')
            mfD.inv{1}.datareg.sensors = grad;
        end
        if isfield(mfD.inv{1},'forward')
            voltype = mfD.inv{1}.forward.voltype;
            mfD.inv{1}.forward = [];
            mfD.inv{1}.forward.voltype = voltype;
            mfD = spm_eeg_inv_forward(mfD,1);
        end
    end
    mfD.save();
end

%-Complete
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));


