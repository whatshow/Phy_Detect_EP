classdef ConstellationEstimIn_Real 
    % ConstellationEstimIn:  Circular Constellation scalar input estimation function
    
    properties
        % Prior mean and variance
        mean0;  % Mean
        var0;   % Variance
        M;      % M-QAM
        x0;     % all constellation points
        % True indicates to compute output for max-sum
        maxSumVal = false;
    end
    
    methods
        % Constructor
        function obj = ConstellationEstimIn_Real(modType, mod_size, mean0, var0, maxSumVal)
            % obj = obj@EstimIn;
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj.M = 2^mod_size ;
                obj.mean0 = mean0 ;
                obj.var0 = var0 ;
                obj.maxSumVal = 0 ;
                
                switch modType
                    case 'QAM',
                        normal_scal = 1/sqrt((2/3)*(2^mod_size-1)) ;   % normalization scale: QPSK = 1/sqrt(2) 16QAM = 1/sqrt(10)
%                         Mod = modem.pammod('M', (2^(mod_size/2)));
                        obj.x0 = normal_scal*pammod([0:(2^(mod_size/2))-1],(2^(mod_size/2)));
                    case 'PSK',
                        normal_scal = 1 ;
                        Mod = modem.pskmod(2^mod_size);
                end
                 
%                 Mod.symbolorder = 'gray';
%                 obj.x0 = normal_scal * Mod.Constellation;

                if (nargin >= 5)
                    if (~isempty(maxSumVal))
                        obj.maxSumVal = maxSumVal;
                    end
                end
            end
        end
        
        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(obj)
            mean0 = obj.mean0;
            var0  = obj.var0;
            obj.maxSumVal = 0 ;
            valInit = 0;
        end
        
        % Circular AWGN estimation function
        % Provides the mean and variance of a variable u
        % from an observation v = u + w, w = CN(0,wvar)
        %
        % function [umean, uvar, val] = estim(obj, v, wvar)
        function [umean, uvar,pxrSum,AllSymbol,pxr] = estim(obj, v, wvar)
            % Get prior
            %umean0 = obj.mean0;
	        %uvar0 = max(eps,obj.var0); % avoid zero variances! 
            %M = obj.M ;
            
            % Compute posterior mean and variance
            AllSymbol = obj.x0 ;
            [size_n,size_m] = size(v) ; 
            v_temp = reshape(v, size_m*size_n, 1) ;
            wvar_temp = reshape(wvar, size_m*size_n, 1) ;
                        
            % logpxr = bsxfun(@times, -1./(2*wvar), abs(bsxfun(@minus, v, AllSymbol)).^2) ;
            % logpxr = bsxfun(@minus, logpxr, max(logpxr, [], 2) );
            % pxr = exp(logpxr);
            % pxr = bsxfun(@rdivide, pxr, sum(pxr,2) ) ; 
            % umean = sum(bsxfun(@times, pxr, AllSymbol), 2) ; 
            % uvar = sum(pxr .* abs(bsxfun(@minus, umean, AllSymbol)).^2, 2) ;
            
            % The above expressions are easily to understand 
            % speed up the peformance ---- 
            if (~obj.maxSumVal)
                logpxr = bsxfun(@times, -1./(2*wvar_temp), abs(bsxfun(@minus, v_temp, AllSymbol)).^2) ;
                logpxr = bsxfun(@minus, logpxr, max(logpxr, [], 2) );
                pxr = exp(logpxr);
                pxrSum = max(sum(pxr,2), eps) ;
                umean_temp = sum(bsxfun(@times, pxr, AllSymbol), 2)./pxrSum ;
                uvar_temp = sum(pxr .* abs(bsxfun(@minus, umean_temp, AllSymbol)).^2, 2)./pxrSum ;
            else
                logpxr = bsxfun(@times, -1./(2*wvar_temp), abs(bsxfun(@minus, v_temp, AllSymbol)).^2) ;
                [xxx, demoIdx] = max(logpxr, [], 2) ;
                umean_temp = AllSymbol(demoIdx) ;
                uvar_temp = wvar;
            end
           
            umean = reshape(umean_temp, size_n, size_m) ;
            uvar = reshape(uvar_temp, size_n, size_m) ;  
            % val = reshape(val, size_n, size_m) ;  
        end
                
    end
    
end

