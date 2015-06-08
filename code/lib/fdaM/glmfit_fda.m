function [bb,dev,stats] = glmfit_fda(x, y, distr, wt, covariate, dataind)
%GLM_FDA Fits a generalized linear model with regularization. 
%  Arguments:
%
%      'covariate' - a vector to use as an additional predictor variable, but
%         with a coefficient value fixed at 1.0.
%
%      'weights' - a vector of prior weights, such as the inverses of the
%         relative variance of each observation.
%
%      'dataind' - a vector of 1's and 0's to indicate which rows of 
%         X and y correspond to actual data, as opposed to regularization 
%         observations.  The default if it is empty is all 1's.  The
%         weights for these values are iteratively modified, while those
%         for the remainder are kept at 1.0.
%
%   Returns:
%       'dfe'       degrees of freedom for error
%       's'         theoretical or estimated dispersion parameter
%       'sfit'      estimated dispersion parameter
%       'se'        standard errors of coefficient estimates B
%       'coeffcorr' correlation matrix for B
%       'covb'      estimated covariance matrix for B
%       't'         t statistics for B
%       'p'         p-values for B
%       'resid'     residuals
%       'residp'    Pearson residuals
%       'residd'    deviance residuals
%       'resida'    Anscombe residuals

%   Last modified 1 February 2013 by Jim Ramsay

%--------------------------------------------------------------------------
%                    Check arguments
%--------------------------------------------------------------------------

if nargin < 3
    error('Number of arguments is less than 3.');
end

%  set default argument values

if nargin < 4, wt         = []; end
if nargin < 5, covariate  = []; end
if nargin < 6, dataind    = []; end
    
N = size(X,1);

%--------------------------------------------------------------------------
%                Set distribution-specific defaults.
%--------------------------------------------------------------------------

switch distr
%--------------------------------------------------------------------------
case 'normal'
    stdFn = @(mu)   ones(size(mu));
    devFn = @(mu,y) (y - mu).^2;
    link  = @(mu)   mu;
    dlink = @(mu)   ones(size(mu));
    ilink = @(eta)  eta;
%--------------------------------------------------------------------------
case 'binomial'
    if size(y,2) == 1
        % N will get set to 1 below
        if any(y < 0 | y > 1)
            error(['For the binomial distribution, ', ...
                   'Y must be a binary vector or\na matrix with ', ...
                   'two columns with the number of trials in the ', ...
                   'second column.']);
        end
    elseif size(y,2) == 2
        y(y(:,2)==0,2) = NaN;
        N = y(:,2);
        y = y(:,1) ./ N;
        if any(y < 0 | y > 1)
            error(['Y must contain values in the interval [0,N] ', ...
                   'for the binomial distribution.']);
        end
    else
        error(['Y must be a two column matrix or a vector for ', ...
               'the binomial distribution.']);
    end
    stdFn = @(mu)   sqrt(mu).*sqrt(1-mu) ./ sqrt(N);
    devFn = @(mu,y) 2*N.*(y.*log((y+(y==0))./mu) + ...
                    (1-y).*log((1-y+(y==1))./(1-mu)));
    link  = @(mu)   log(mu ./ (1-mu));
    dlink = @(mu)   1 ./ (mu .* (1-mu));
    loBnd = log(eps); 
    upBnd = -loBnd;
    ilink = @(eta)  1 ./ (1 + exp(-constrain(eta,loBnd,upBnd)));
%--------------------------------------------------------------------------
case 'poisson'
    if any(y < 0)
        error(['Y must contain non-negative values for ', ...
            'the Poisson distribution.']);
    end
    stdFn = @(mu)   sqrt(mu);
    devFn = @(mu,y) 2*(y .* (log((y+(y==0)) ./ mu)) - (y - mu));
    link  = @(mu)   log(mu);
    dlink = @(mu)   1 ./ mu;
    lowerBnd = log(eps); upperBnd = -lowerBnd;
    ilink = @(eta)  exp(constrain(eta,lowerBnd,upperBnd));
%--------------------------------------------------------------------------
case 'gamma'
    if any(y <= 0)
        error('glm_fit:BadData', ...
              ['Y must contain positive values for ', ...
               'the gamma distribution.']);
    end
    % keep mu = ilink(eta) in [eps, 1/eps];
    stdFn = @(mu) mu;
    devFn = @(mu,y) 2*(-log(y ./ mu) + (y - mu) ./ mu);
    link  = @(mu)  1 ./ mu;
    dlink = @(mu) -1 ./ mu.^2;
    lowerBnd = eps; upperBnd = 1/lowerBnd;
    ilink = @(eta) 1 ./ constrain(eta,lowerBnd,upperBnd);
%--------------------------------------------------------------------------
case 'inverse gaussian'
    if any(y <= 0)
        error('glm_fit:BadData', ...
              ['Y must contain positive values for ', ...
               'the inverse Gaussian distribution.']);
    end
    % linkExponent==0 (equivalent to 'log') has been weeded out already
    % keep eta = link(mu) in [eps, 1/eps];
    % keep mu = ilink(eta) in [eps, 1/eps];
    stdFn = @(mu) mu.^(3/2);
    devFn = @(mu,y) ((y - mu)./mu).^2 ./  y;
    Expon = -2;
    lowerBnd = eps.^min(abs(linkExponent),1); upperBnd = 1/lowerBnd;
    link  = @(mu)  constrain(mu,lowerBnd,upperBnd).^Expon;
    dlink = @(mu)  Expon * mu.^(Expon-1);
    lowerBnd = eps.^min(abs(1/linkExponent),1); upperBnd = 1/lowerBnd;
    ilink = @(eta) constrain(eta,lowerBnd,upperBnd).^(1/Expon);
%--------------------------------------------------------------------------
otherwise
    error('Distribution name is invalid.');
end

%  check DATAIND to set to indicate all rows correspond to actual data

if isempty(dataind)
    dataind = ones(n,1);
else
    if length(dataind) > n
        error('DATAIND is too long.');
    end
end
dataind   = (dataind == 1);

[n,ncolx] = size(x);

%  define default weight vector PWTS and check for positivity

if isempty(wt)
    wt = 1;
elseif any(wt == 0)
    n = n - sum(wt == 0);
end

%--------------------------------------------------------------------------
%                    Set up for iterations
%--------------------------------------------------------------------------

% Initialize mu and eta from y.

mu  = startingVals(distr,y,N);
mu(~dataind)  = y(~dataind);
eta = linkFun(mu);
eta(~dataind) = y(~dataind);

iter     = 0;
iterLim  = 100;
warned   = false;
seps     = sqrt(eps);
convcrit = 1e-6;
b        = zeros(p,1,dataClass);

% Enforce limits on mu to guard against an inverse link that doesn't map 
% into the support of the distribution.

switch distr
case 'binomial'
    % mu is a probability, so order one is the natural scale, and eps is a
    % reasonable lower limit on that scale (plus it's symmetric).
    muLims = [eps(dataClass) 1-eps(dataClass)];
case {'poisson' 'gamma' 'inverse gaussian'}
    % Here we don't know the natural scale for mu, so make the lower limit
    % small.  This choice keeps mu^4 from underflowing.  No upper limit.
    muLims = realmin(dataClass).^.25;
end

%--------------------------------------------------------------------------
%                    start of GLM iteration loop
%--------------------------------------------------------------------------

sqrtwt = sqrt(wt);

while iter <= iterLim
    iter = iter+1;

    % Compute adjusted dependent variable for least squares fit
    
    Deta = dlinkFun(mu);
    Deta(~dataind) = 1;
    z    = eta + (y - mu) .* Deta;

    % Compute IRLS weights the inverse of the variance function
    
    sqrtw = sqrtwt ./ (abs(Deta) .* stdFn(mu));
    
    % Replace weights for non-data rows by 1's
    
    sqrtw(~dataind) = 1;

    % Check sqrtw
    
    wtol = max(sqrtw)*eps(dataClass)^(2/3);
    t = (sqrtw < wtol);
    if any(t)
        t = t & (sqrtw ~= 0);
        if any(t)
            sqrtw(t) = wtol;
            if ~warned
                warning(...
                  ['Weights are ill-conditioned.   Data may be badly ' ...
                   'scaled, or\nthe link function may be inappropriate.']);
            end
            warned = true;
        end
    end

    % Compute coefficient estimates for this iteration - the IRLS step
    
    b_old = b;
    ytmp  = z - covariate;
    yw    = ytmp .* sqrtw;
    xw    = x .* sqrtw(:,ones(1,p));
    [Q,R] = qr(xw,0);
    b     = R \ (Q'*yw);
    eta   = covariate + x*b;
    mu    = ilinkFun(eta);
    mu(~dataind) = eta(~dataind);

    % Force mean in bounds, in case the link function is a wacky choice
    
    switch distr
    case 'binomial'
        if any(mu(dataind) < muLims(1) | muLims(2) < mu(dataind))
        mu(dataind) = max(min(mu(dataind),muLims(2)),muLims(1));
        end
    case {'poisson' 'gamma' 'inverse gaussian'}
        if any(mu(dataind) < muLims(1))
        mu(dataind) = max(mu(dataind),muLims(1));
        end
    end

    % Check stopping conditions
    
    if (~any(abs(b-b_old) > convcrit * max(seps, abs(b_old)))), break; end
end

%--------------------------------------------------------------------------
%                    end of GLM iteration loop
%--------------------------------------------------------------------------

if iter > iterLim
    warning('Iteration limit reached.');
end

bb = zeros(ncolx,1,dataClass); bb(perm) = b;

if nargout > 1
    % Sum components of deviance to get the total deviance.
    di  = devFn(mu,y);
    di(~dataind) = (y(~dataind) - mu(~dataind)).^2;
    dev = sum(wt .* di);
end

%--------------------------------------------------------------------------
%                Compute additional statistics 
%--------------------------------------------------------------------------

if nargout > 2
    % Compute the sum of squares used to estimate dispersion, and the
    % Anscombe residuals.
    switch(distr)
    case 'normal'
        ssr = sum(wt(dataind) .* (y(dataind) - mu(dataind)).^2);
        anscresid = y(dataind) - mu(dataind);
    case 'binomial'
        ssr = sum(wt(dataind) .* (y(dataind) - mu(dataind)).^2 ./ ...
            (mu(dataind) .* (1 - mu(dataind)) ./ N(dataind)));
        t = 2/3;
        anscresid = beta(t,t) * ...
            (betainc(y(dataind),t,t)-betainc(mu(dataind),t,t)) ./ ...
            ((mu(dataind).*(1-mu(dataind))).^(1/6) ./ sqrt(N(dataind)));
    case 'poisson'
        ssr = sum(wt(dataind) .* (y(dataind) - mu(dataind)).^2 ./ ...
            mu(dataind));
        anscresid = 1.5 * ((y(dataind).^(2/3) - mu(dataind).^(2/3)) ./ ...
            mu(dataind).^(1/6));
    case 'gamma'
        ssr = sum(wt .* ((y - mu) ./ mu).^2);
        anscresid = 3 * (y.^(1/3) - mu.^(1/3)) ./ mu.^(1/3);
    case 'inverse gaussian'
        ssr = sum(wt(dataind) .* ((y(dataind) - mu(dataind)) ./ ...
            mu(dataind).^(3/2)).^2);
        anscresid = (log(y(dataind)) - log(mu(dataind))) ./ mu(dataind);
    end

    % Compute residuals, using original count scale for binomial
    
    if (isequal(distr, 'binomial'))
        resid = (y(dataind) - mu(dataind)) .* N(dataind);
    else
        resid  = y(dataind) - mu(dataind);
    end

    dfe = max(n - p, 0);
    stats.beta = bb;
    stats.dfe  = dfe;
    if dfe > 0
        stats.sfit = sqrt(ssr / dfe);
    else
        stats.sfit = NaN;
    end
    if isequal(estdisp,'off')
        stats.s = 1;
        stats.estdisp = false;
    else
        stats.s = stats.sfit;
        stats.estdisp = true;
    end

    % Find coefficient standard errors and correlations
    if ~isnan(stats.s) % dfe > 0 or estdisp == 'off'
        RI = R\eye(p);
        C = RI * RI';
        if isequal(estdisp,'on'), C = C * stats.s^2; end
        se = sqrt(diag(C)); se = se(:);   % insure vector even if empty
        stats.covb = zeros(ncolx,ncolx,dataClass);
        stats.covb(perm,perm) = C;
        C = C ./ (se * se');
        stats.se = zeros(ncolx,1,dataClass); stats.se(perm) = se;
        stats.coeffcorr = zeros(ncolx,ncolx,dataClass);
        stats.coeffcorr(perm,perm) = C;
        stats.t = NaN(ncolx,1,dataClass); stats.t(perm) = b ./ se;
        if isequal(estdisp,'on')
            stats.p = 2 * tcdf(-abs(stats.t), dfe);
        else
            stats.p = 2 * normcdf(-abs(stats.t));
        end
    else
        stats.se = NaN(size(bb),class(bb));
        stats.coeffcorr = NaN(length(bb),class(bb));
        stats.t = NaN(size(bb),class(bb));
        stats.p = NaN(size(bb),class(bb));
    end

    stats.resid  = statinsertnan(wasnan, resid);
    stats.residp = ...
            statinsertnan(wasnan, (y - mu) ./ (stdFn(mu) + (y==mu)));
    stats.residd = statinsertnan(wasnan, sign(y - mu) .* sqrt(max(0,di)));
    stats.resida = statinsertnan(wasnan, anscresid);
end



