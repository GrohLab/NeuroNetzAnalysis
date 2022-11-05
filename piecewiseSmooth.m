function [pw_yhat] = piecewiseSmooth(signal, scale, m)

pw_yhat = nan(length(scale), length(signal));
for s = 1:length(scale)
    idxs = 1:scale(s);
    w = linspace(0, 1, round(scale(s)/2));
    while max(idxs) <= length(signal)
        [polft, ~]=detrend_profile(m,idxs,signal(idxs));
        pw_yhat(s,idxs) = sum([...
            pw_yhat(s,idxs).*[(1-w),nan(1,round(scale(s)/2))];...
            polft.*[w,ones(1,round(scale(s)/2))]],1,'omitnan');
        idxs = idxs + round(scale(s)/2);
    end
end
end

