function kpi = computeKeypoints(V, F, curve, scale)
    n = 10;
    minLength = 0.03 / scale;   % 3 cm
    maxLength = 0.15 / scale;   % 15 cm
    

    [gk, ~, path] = geodesicCurvature(V, F, curve);
    P = [vertcat(cell2mat(path).x), vertcat(cell2mat(path).y) vertcat(cell2mat(path).z)];
    
    edge = diff(P, [], 1);
    eLen = vecnorm(edge, 2, 2);
    L = sum(eLen);
    cumLen = [0; cumsum(eLen)];
    minLen = max(L/n, minLength);
    badLenIdx = diff(cumLen)<1e-9;
    gk(badLenIdx) = [];
    cumLen(badLenIdx) = [];
    
    signal = abs(movmean(gk, L/(2*n), 'SamplePoints', cumLen));
    
    [peakVal, candidates] = findpeaks(signal, 'MinPeakHeight', mean(signal));
    
    kpi = [1; numel(gk)];
    
    
    % sort candidates by the peak values
    [~, idx] = sort(peakVal, 'descend');
    candidates = candidates(idx);
    
    % greedily choose candidate points
    for i=1:numel(candidates)
        curCand = candidates(i);
        distCand = abs(cumLen(curCand) - cumLen(kpi));
        if min(distCand) > minLen
            kpi = [kpi; curCand];
        end
    end
    

    kpi = sort(kpi);
    while 1
        lenDiff = diff(cumLen(kpi));
        longSections = find(lenDiff > maxLength);
        if numel(longSections)==0
            break;
        end
        
        newPtIdx = [];
        for i=1:numel(longSections)
            sectionLen = lenDiff(longSections(i));
            numNewPts = floor(sectionLen / maxLength);
            newPtLocs = cumLen(kpi(longSections(i))) + (1:numNewPts)/(numNewPts+1) * sectionLen;
            n = numel(newPtLocs);
            for j=1:n
                newPtIdx(end+1) = find(cumLen > newPtLocs(j), 1);
            end
        end
        kpi = sort([kpi; newPtIdx']);
    end

    pts = getPoints(V, F, curve);
    
    [~, kpi] = pdist2(pts, P(kpi, :), 'squaredeuclidean', 'Smallest', 1);
    kpi = sort(kpi)';
end
