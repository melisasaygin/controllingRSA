function [repairddata, CompleetTime, repairdIndex] = repair_holes(t, data)
    % Sort the time vector if not already sorted
    if ~issorted(t)
        [ttsort, index] = sort(t);
        repairddata = data(index);
    else
        ttsort = t;
        repairddata = data;
    end

    % Adjust time to start from 0
    starttimeTemp = ttsort(1) - 1;
    ttsort = ttsort - starttimeTemp;

    % Remove duplicates
    [ttsort, ~, ic] = unique(ttsort, 'stable');
    repairddata = accumarray(ic, repairddata, [], @mean);

    % Fill holes using interpolation
    CompleetTime = ttsort(1):ttsort(end);
    repairddata = interp1(ttsort, repairddata, CompleetTime, 'linear', 'extrap');

    % Identify missing time indices and adjust to original scale
    missingTimes = setdiff(CompleetTime, ttsort) + starttimeTemp;
    repairdIndex = missingTimes;
    CompleetTime = CompleetTime + starttimeTemp;
end
