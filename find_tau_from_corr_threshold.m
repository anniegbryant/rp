%% --- Helper function: find tau ---
function tau = find_tau_from_corr_threshold(x, m, corr_type, threshold)

    if nargin < 4
        threshold = 0;
    end

    N = length(x);
    max_lag = round(0.2 * N); % limit lag search to 20% of series length

    if strcmpi(corr_type, 'ACF')
        mid     = length(x);
        acfvals = xcorr(x-mean(x),'coeff');
        acf_pos = acfvals(mid+1:mid+max_lag); % positive lags only
        tau = find(acf_pos <= threshold, 1 ,'first');
    else
        % requires hctsa's CO_FirstMin.m
        % Compute first minimum of the (auto)mutual information with Kraskov estimator (type 2)
        tau = CO_FirstMin(x, 'mi', 'kraskov2');
    end
    
    if isempty(tau)
        tau = 1; % fallback
    end
    
    % Cap tau at the max_lag if necessary
    if tau > max_lag
        warning('Tau=%d too large, capping to %d', tau, max_lag);
        tau = max(1,max_lag);
    else
        fprintf('Using Tau=%d\n', tau);
    end
end