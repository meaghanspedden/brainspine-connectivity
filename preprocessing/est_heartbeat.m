function [heartep, beatlen] = est_heartbeat(D, heartind, megind, heartidx, flag, beatlen)

if nargin < 6
    beatlen = [];
end

data = squeeze(D(heartind(heartidx), :));
fs = D.fsample;

heartest = data(1, :);

if flag
    windowSize = fs * 30; % 30 sec
    stepSize = fs * 10; % 10 sec
    numWindows = floor((length(heartest) - windowSize) / stepSize) + 1;
    qualityScores = zeros(1, numWindows);

    for i = 1:numWindows
        idxStart = (i-1)*stepSize + 1;
        idxEnd = idxStart + windowSize - 1;
        segment = heartest(idxStart:idxEnd);
        qualityScores(i) = 1 / var(segment);
    end

    [~, bestWinIdx] = max(qualityScores);
    bestSegmentIdxStart = (bestWinIdx - 1) * stepSize + 1;
    bestSegment = heartest(bestSegmentIdxStart : bestSegmentIdxStart + windowSize - 1);

    heartest = bestSegment;  % Use best segment
    offset = bestSegmentIdxStart - 1;  % Track offset for indexing in D
else
    offset = 0;
end

[b, a] = butter(2, [0.5 40] / (fs / 2), 'bandpass');
heartest = filtfilt(b, a, heartest);

minPeakDist = round(0.6 * fs);

[peakPos, locPos] = findpeaks(heartest, 'MinPeakDistance', minPeakDist);
[peakNeg, locNeg] = findpeaks(-heartest, 'MinPeakDistance', minPeakDist);
peakNeg = -peakNeg;

N = min(50, min(length(peakPos), length(peakNeg)));
topPos = sort(peakPos, 'descend');
topNeg = sort(abs(peakNeg), 'descend');

avgPos = median(topPos(1:N));
avgNeg = median(topNeg(1:N));

if avgNeg > avgPos
    heartest = -heartest;
end

heightthresh = 0.7 * max(abs([avgPos avgNeg]));

[pks, htrigs] = findpeaks(heartest, 'MinPeakDistance', minPeakDist, 'MinPeakHeight', heightthresh);

zScores = zscore(pks);
inliers = abs(zScores) < 2;
htrigs = htrigs(inliers);

segmentDuration = length(heartest) / fs;
hbfreq = length(htrigs) / segmentDuration * 60;
fprintf('Estimated heart rate: %g bpm\n', hbfreq);

if isempty(beatlen)
    beatlen = 2 * floor(median(diff(htrigs)) / 2); % even number
end

heartep = [];

for f = 1:length(htrigs)
    idx = htrigs(f) + offset;
    startIdx = idx - beatlen / 2;
    endIdx = idx + beatlen / 2 - 1;
    if startIdx > 0 && endIdx <= size(D, 2)
        heartep(:, :, f) = D(megind, startIdx:endIdx, 1);
    end
end

heartep = mean(heartep, 3);

figure; plot(heartep');
legend('Estimate of heartbeat (all chans)');

end
