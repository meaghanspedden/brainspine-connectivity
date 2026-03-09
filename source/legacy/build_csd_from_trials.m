function C = build_csd_from_trials(freqdat_tr, trialIdx)
    nCh = numel(freqdat_tr.label);
    C = complex(zeros(nCh,nCh));

    for t = trialIdx(:)'
        for k = 1:size(freqdat_tr.labelcmb,1)
            ch1 = find(strcmp(freqdat_tr.label, freqdat_tr.labelcmb{k,1}));
            ch2 = find(strcmp(freqdat_tr.label, freqdat_tr.labelcmb{k,2}));
            val = freqdat_tr.crsspctrm(t,k);
            C(ch1,ch2) = C(ch1,ch2) + val;
            C(ch2,ch1) = C(ch2,ch1) + conj(val);
        end
    end

    C = C / numel(trialIdx);  % average
end