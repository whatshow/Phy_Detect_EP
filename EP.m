classdef EP < handle
    properties
        constellation {mustBeNumeric}
        constellation_len {mustBeNumeric}
        es = 1;                                 % constellation average power
        beta = 0.2;                             % the percentage for taking values from the previous iteration 
        epsilon {mustBeNumeric} = 2*5e-7;       % the default minimal variance is 2*5e-7 (5e-7 for real only)
        l = 10;                                 % maximal iteration
        early_stop = false;
        early_stop_min_diff = 2e-4;             % the elemental minimal difference for mean and variance in Bayesian estimation (1e-4 for real only)
    end
    methods
        % init
        % @constellation:           the constellation, a vector
        % @beta:                    the percentage for taking values from the previous iteration
        % @epsilon:                 the default minimal variance
        % @l:                       maximal iteration
        % @early_stop:              whether stop early
        % @early_stop_min_diff:     the elemental minimal difference in mean and variance
        function self = EP(constellation, varargin)
            % inputs
            % inputs - constellation
            if ~isvector(constellation)
                error("The constellation must be a vector.");
            else
                % constellation must be a row vector or an 1D vector
                constellation = constellation(:);
                constellation = constellation.';
                self.constellation = constellation;
                self.constellation_len = length(constellation);
            end
            % inputs - constellation average power
            self.es = sum(self.constellation.^2)/self.constellation_len;
            % optional inputs - register 
            inPar = inputParser;
            addParameter(inPar,"beta", self.beta, @(x) isscalar(x)&isnumeric(x));
            addParameter(inPar,"epsilon", self.epsilon, @(x) isscalar(x)&isnumeric(x));
            addParameter(inPar,"l", self.l, @(x) isscalar(x)&isnumeric(x));
            addParameter(inPar,"early_stop", self.early_stop, @(x) isscalar(x)&islogical(x));
            addParameter(inPar,"early_stop_min_diff", self.early_stop_min_diff, @(x) isscalar(x)&isnumeric(x));
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            % optional inputs - take
            self.beta = inPar.Results.beta;
            self.epsilon = inPar.Results.epsilon;
            self.l = inPar.Results.l;
            self.early_stop = inPar.Results.early_stop;
            self.early_stop_min_diff = inPar.Results.early_stop_min_diff;
        end
        
        % detect
        % @y:           the received signal
        % @H:           the channel matrix
        % @No:          the noise (linear) power
        % @sym_map:     false by default. If true, the output will be mapped to the constellation
        function x = detect(self, y, H, No, varargin)
            % inputs
            if isscalar(y) 
                error("The received signal must be a vector.")
            elseif ~isvector(y)
                error("The received signal must be a vector.")
            end
            if isscalar(H) 
                error("The channel must be a matrix.")
            elseif ~ismatrix(H)
                error("The channel must be a matrix.")
            end
            [y_num, x_num] = size(H);
            if y_num ~= length(y)
                error("The channel row number does not equal to the signal number.");
            end
            if y_num < x_num
                error("The channel is a correlated channel.")
            end
            if ~isscalar(No)
                error("The noise power must be a scalar.");
            end
            % optional inputs - register 
            inPar = inputParser;
            addParameter(inPar,"sym_map", false, @(x) isscalar(x)&islogical(x)); 
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            % optional inputs - take
            sym_map = inPar.Results.sym_map;
            
            % constant values
            Ht = H';
            Hty = Ht*y;
            HtH = Ht*H;
            
            % iterative detection
            gamma = zeros(x_num, 1);
            lamda = (1/self.es)*ones(x_num, 1);
            if self.early_stop
                lmmse_mu_prev = zeros(x_num, 1);
                lmmse_Sigma_prev = zeros(x_num, x_num);
            end
            for iter_id = 1:self.l
                % linear MMSE
                lmmse_Sigma  = inv(HtH + No*diag(lamda));       % covariance matrix
                lmmse_mu = lmmse_Sigma*(Hty + No*gamma);        % mean vector
                
                % early stop
                if self.early_stop
                    if iter_id > 1 && sum(abs(lmmse_mu_prev - lmmse_mu) < self.early_stop_min_diff ) == x_num && sum(abs(lmmse_Sigma_prev - lmmse_Sigma) < self.early_stop_min_diff, "all") == x_num*x_num
                        break;
                    end
                    % early stop - store mean and covariance
                    lmmse_mu_prev = lmmse_mu;
                    lmmse_Sigma_prev = lmmse_Sigma;
                end
                
                % cavity marginal
                sigmai2 = No*diag(real(lmmse_Sigma));
                hi2 = sigmai2./(1-sigmai2.*lamda);
                hi2 = max(hi2, eps);
                ti = hi2.*(lmmse_mu./sigmai2 - gamma);
                
                % Bayessian Estimation
                pxyPdfExpPower = -1./(2*hi2).*abs(repmat(ti, 1, self.constellation_len) - repmat(self.constellation, x_num, 1)).^2;
                pxypdfExpNormPower = pxyPdfExpPower - max(pxyPdfExpPower, [], 2);   % make every row the max power is 0
                pxyPdf = exp(pxypdfExpNormPower);
                % Bayessian Estimation - Calculate the coefficient of every possible x to make the sum of all
                pxyPdfCoeff = 1./sum(pxyPdf, 2);
                pxyPdfCoeff = repmat(pxyPdfCoeff, 1, self.constellation_len);
                % Bayessian Estimation - PDF normalisation
                pxyPdfNorm = pxyPdfCoeff.*pxyPdf;
                % Bayessian Estimation - calculate the mean and variance
                mu_pi = sum(pxyPdfNorm.*self.constellation, 2);
                mu_pi_mat = repmat(mu_pi, 1, self.constellation_len);
                sigma_pi = sum(abs(mu_pi_mat - self.constellation).^2.*pxyPdfNorm, 2);
                sigma_pi = max(sigma_pi, self.epsilon);
                
                % update
                lamda_new = (hi2-sigma_pi)./sigma_pi./hi2;
				gamma_new = (mu_pi.*hi2-sigma_pi.*ti )./sigma_pi./hi2 ;  
                % update - avoid Negative Index
                lamda_new_neg_idx = lamda_new < 0;
                lamda_new(lamda_new_neg_idx) = lamda(lamda_new_neg_idx);
                gamma_new(lamda_new_neg_idx) = gamma(lamda_new_neg_idx);
                % update - convex combination
                lamda = (1-self.beta)*lamda + self.beta*lamda_new;
                gamma = (1-self.beta)*gamma + self.beta*gamma_new;
            end
            x = mu_pi;
            
            % hard estimation
            if sym_map
                x = self.symmap(x);
            end
        end
        
        % symbol mapping (hard)
        % @syms: a vector of symbols
        function syms_mapped = symmap(self, syms)
            if ~isvector(syms)
                error("Symbols must be into a vector form to map.");
            end
            % the input must be a column vector
            is_syms_col = iscolumn(syms);
            syms = syms(:);
            syms_dis = abs(syms - self.constellation).^2;
            [~,syms_dis_min_idx] =  min(syms_dis,[],2);
            syms_mapped = self.constellation(syms_dis_min_idx);
            if is_syms_col
                syms_mapped = syms_mapped(:);
            end
        end
    end
end