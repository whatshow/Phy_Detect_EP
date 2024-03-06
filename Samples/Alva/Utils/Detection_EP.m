%% EP Detection 
% 1.This uses EP algorithm to detect the M-ary modulation symbols (for M-QAM, M-PSK, etc.)
% 2.H must be known and each X is fixed, so H*x is a constant in each symbol period
% -------------------------------------------------------------------------
% OUTPUTS
%
% @syms: all dectected symbols based on the input
% -------------------------------------------------------------------------
% INPUTS
%
% @sympool:         The pool has all possible symbols's values              (This MUST be a row matrix)
% @H:               Channel estimation 
% @y:               Received Signals                                        (This MUST a column matrix)
% @y_var:           Received Signals' variance (This must be a scalar)
% @ep_iter_times:   How many times we should do EP likelihood calculation
% @beta:            Convex combination, usually is 0.2 (if no input, beta will be 0.2)

function [syms] = Detection_EP(sympool, H, y, y_var, ep_iter_times, varargin)
    %% load multiple inputs
    % define the default value
    beta = 0.2;
    pxy_pred_var_min = 1e-13;
    pxy_pred_mean_diff = 1e-4;
    % load multiple inputs
    inPar = inputParser;
    % Register names
    addParameter(inPar,'Beta', beta, @(x)isnumeric(x)&&isscalar(x) );
    addParameter(inPar,'MinVariance', pxy_pred_var_min, @(x)isnumeric(x)&&isscalar(x) );
    % Allow unmatched cases
    inPar.KeepUnmatched = true; 
    % Allow capital or small characters
    inPar.CaseSensitive = false;
    % Try to load those inputs 
    parse(inPar, varargin{:}); 
    if isempty(find(inPar.UsingDefaults == "Beta", 1))                     % if input exists, use it 
        beta = inPar.Results.Beta;
    end
    if isempty(find(inPar.UsingDefaults == "MinVariance", 1))              % if input exists, use it 
        pxy_pred_var_min = inPar.Results.MinVariance;
    end
    
    % -------------------------------------------------------------------------
    % DEV
    % -------------------------------------------------------------------------
    
    % Get Symbolpool capacity
    sympool_len = length(sympool);
    % retrieve tx antenna number and rx antenna number
    [~, rx_num] = size(H);
    
    % -------------------------------------------------------------------------
    % EP algorithm
    % -------------------------------------------------------------------------
    % P(x|y) = P(y|x)P(x)
    % Param Config - Every Iteration
    px_Lambda = 2*ones(rx_num, 1);                                % p(x) variance reciprocal
    px_gamma = zeros(rx_num, 1);                                % p(x) mean is gamma/Lambda
    px_Lambda_new = ones(rx_num, 1);                            % p(x) variance reciprocal (new iteration)
    px_gamma_new = zeros(rx_num, 1);                            % p(x) mean (new iteration)
  
    
    Ht = H';
    HtH = Ht*H;
    
    % Try 'ep_iter_times' times to estimate the P(x|y)
    for i = 1:1:ep_iter_times
        % Calculate P(x|y) current mean & variance
        pxy_covar_TMP = inv(HtH + y_var*diag(px_Lambda));
        pxy_covar = y_var*pxy_covar_TMP;
        pxy_mean = pxy_covar_TMP*(Ht*y + y_var*px_gamma);
        pxy_var = real(diag(pxy_covar));                              % Get Variance from Covariance (The diagonal is the covariance)

        
        % p(y|x) Cavity Distribution - Mean & Variance: p(y|x) = p(x|y)/p(x) 
        pyx_var = max(pxy_var./(1 - pxy_var.* px_Lambda), eps);
        pyx_mean = pyx_var.*(pxy_mean./pxy_var - px_gamma);
        
        % Find the mean of Prediction P(x|y) by P(y|x)*P(x), P(x) is uniform distributed
        pyx_mean_mat = repmat(pyx_mean, 1, sympool_len);                                    % copy 'sympool_len' colums of 'pyx_mean'
        sympool_mat = repmat(sympool, rx_num, 1);                                           % sympool for all Tx antennas
        % Calculate P(x|y) for every x (For complex values, we use absolute diference)
        pxy_pred_pdf_power = -1./(2*pyx_var).*abs(pyx_mean_mat - sympool_mat).^2;               % e^pxy_pred_pdf_power (for all Tx and all possible symbols)
        pxy_pred_pdf_power = pxy_pred_pdf_power - max(pxy_pred_pdf_power, [], 2);           % subtract the maximum pdf_power on every possible selection of each row (for all Tx and all possible symbols)
        pxy_pred_pdf = exp(pxy_pred_pdf_power);                                             % PDF (for all Tx and all possible symbols)
        pxy_pred_pdf_sum = sum(pxy_pred_pdf, 2);
        pxy_pred_pdf_coeff = 1./pxy_pred_pdf_sum;
        pxy_pred_pdf_coeff_mat = repmat(pxy_pred_pdf_coeff, 1, sympool_len);                % PDF coeff(a column vector) for every symbols
        pxy_pred_pdf_norm = pxy_pred_pdf.*pxy_pred_pdf_coeff_mat;                           % Norm PDF (The sum of each row is 1)
        % If the difference between pxy_pred_mean old and new value are less than 10^-4
%         if i > 1
%             if abs(pxy_pred_mean - sum(pxy_pred_pdf_norm.*sympool_mat, 2)) < pxy_pred_mean_diff
%                 break;
%             end
%         end
        pxy_pred_mean = sum(pxy_pred_pdf_norm.*sympool_mat, 2);                             % ONLY for all Tx
        pxy_pred_mean_mat = repmat(pxy_pred_mean, 1, sympool_len);                          % mean (for all Tx and all possible symbols)
        pxy_pred_var = sum(abs(sympool_mat - pxy_pred_mean_mat).^2.*pxy_pred_pdf_norm, 2);  % ONLY for all Tx
        % P(x|y) prediction variance at least should be 10^-13
        pxy_pred_var = max(pxy_pred_var_min, pxy_pred_var);
        
        
        % Create new Lambda and Gamma
        px_Lambda_new_origin = (pyx_var - pxy_pred_var)./pxy_pred_var./pyx_var;
        px_gamma_new_origin  = (pxy_pred_mean.*pyx_var - pyx_mean.*pxy_pred_var)./pyx_var./pxy_pred_var;
        
        % Replace negative parts in new Lambda and Gamma with old values        
        px_Lambda_new_nn = px_Lambda_new_origin;
        px_gamma_new_nn = px_gamma_new_origin;
        % Find negative Idx
        negIdx = px_Lambda_new_nn < 0;
        px_Lambda_new_nn(negIdx) = px_Lambda(negIdx);
        px_gamma_new_nn(negIdx) = px_gamma(negIdx);
        
        % Smooth new Lambda and Gamma
        px_Lambda_new = beta*px_Lambda + (1 - beta)*px_Lambda_new_nn;
        px_gamma_new = beta*px_gamma + (1 - beta)*px_gamma_new_nn;
        
        % Use new Lambda and Gamma
        px_Lambda = px_Lambda_new;
        px_gamma = px_gamma_new;
    end
    % Give the guess symbols from Predication P(x|y)
    syms = pxy_pred_mean;
    %syms = pyx_mean;
end