function y = embed(varargin)
%EMBED   Creates embedding vector using time delay embedding
%     Y=EMBED(X,M,T) creates the embedding vector Y from the time
%     series X using a time delay embedding with dimension M and
%     delay T. The resulting embedding vector has length N-T*(M-1),
%     where N is the length of the original time series.
%
%     Y=EMBED(X,M,'acf') automatically selects tau as the first
%     zero crossing of the autocorrelation function.
%
%     Y=EMBED(X,M,'mi') automatically selects tau as the first
%     minimum of mutual information (requires CO_FirstMin).
%
%     Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         y = embed(x,2,'acf'); % embed with auto-selected tau
%         plot(y(:,1),y(:,2))

%% check input and output arguments
narginchk(1,3)
nargoutchk(0,1)

%% defaults
tau = 1; 
m   = 1;

%% parse arguments
if nargin > 2
    tau = varargin{3};
end
if nargin > 1
    m = varargin{2};
end    
x = varargin{1}; % input vector

%% auto-tau case
if ischar(tau)
    switch lower(tau)
        case 'acf'
            tau = find_tau_from_corr_threshold(x, m, 'ACF', 0);
        case 'mi'
            tau = find_tau_from_corr_threshold(x, m, 'MI');
        otherwise
            error('Unknown tau option: use numeric value, ''acf'', or ''mi''.')
    end
end

%% length checks
N = size(x,1) - (m-1)*tau; 
if N <= 1 % check if row vector (transform it to column vector)
   x = x(:);
   N = length(x) - (m-1)*tau; 
end

d = size(x,2); 
if d > N
   error('Number of columns should be smaller than the length of the time series.')
end

%% create embedding vector 
if size(x,2) == 1 % input vector is one-column time series
    y = buffer(x,N,N-tau,'nodelay');
else % input vector is multi-column time series
    y = zeros(N, d*m);
    for i = 1:d
       y(:, (i-1)*m+(1:m) ) = buffer(x(:,i),N,N-tau,'nodelay');
    end
end
end % end of main embed function

