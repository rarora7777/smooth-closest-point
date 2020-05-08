function kpi = computeKeypoints(V, F, curve, scale)
    n = 10;
    minLength = 0.03 / scale;   % 3 cm
    
%     fIdx = curve.FI;
%     bary = curve.B;
%     bary(bary < 0) = bary;
%     bary = bary ./ sum(bary, 2);
    
    [gk, ~, path] = geodesicCurvature(V, F, curve);
    P = [vertcat(cell2mat(path).x), vertcat(cell2mat(path).y) vertcat(cell2mat(path).z)];
    
    edge = diff(P, [], 1);
    eLen = vecnorm(edge, 2, 2);
%     vLen = ([0; eLen] + [eLen; 0])/2;
    L = sum(eLen);
    cumLen = [0; cumsum(eLen)];
    minLen = max(L/n, minLength);
    badLenIdx = diff(cumLen)<1e-9;
    gk(badLenIdx) = [];
    cumLen(badLenIdx) = [];
    
    signal = abs(movmean(gk, L/(2*n), 'SamplePoints', cumLen));
    
    [peakVal, candidates] = findpeaks(signal, 'MinPeakHeight', mean(signal));
    
    midPt = 1+find(cumLen > L/2, 1);
    kpi = [1; midPt; numel(gk)];
    
    peakVal(candidates <= midPt) = [];
    candidates(candidates <= midPt) = [];
    
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
    
%     prevPt = midPt;
%     kpi = midPt;
%     
%     while cumLen(prevPt) < L*(n-1)/n
%         cdIdx = find(cumLen(candidates) > cumLen(prevPt) + L/n, 1);
%         if isempty(cdIdx)
%             break;
%         end
%         kpi = [kpi; candidates(cdIdx)];
%         prevPt = kpi(end);
%     end
    
%     kpi = [1; kpi; numel(gk)];
    
    pts = getPoints(V, F, curve);
    
    [~, kpi] = pdist2(pts, P(kpi, :), 'squaredeuclidean', 'Smallest', 1);
    kpi = sort(kpi)';
end
