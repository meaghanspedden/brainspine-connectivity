function C = csd_from_trialset(freqdat_tr, trialset)
    labels = freqdat_tr.label;
    labelcmb = freqdat_tr.labelcmb;

    nCh = numel(labels);
    C = complex(zeros(nCh,nCh));

    % average cross-spectra across trials
    csd_avg = mean(freqdat_tr.crsspctrm(trialset, :), 1);  % 1 x nPairs

    for k = 1:size(labelcmb,1)
        ch1 = find(strcmp(labels, labelcmb{k,1}));
        ch2 = find(strcmp(labels, labelcmb{k,2}));
        val = csd_avg(k);
        C(ch1,ch2) = val;
        C(ch2,ch1) = conj(val);
    end
end