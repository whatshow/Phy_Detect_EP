%% EP detection 
% <Author>
% Alva Kosasih
% <INPUTS>
% @y: a column vector for received signal
% @H: the channel matrix
% @noiseLevel: the noise electrical level in voltage (not in dB). noiseLevel = 10^(-SNR/10);
% @modType: 'PSK' or 'QAM'
% @mod_size: the number bits for 1 constellation point in any M-ary Modulation, which is log2(M). For example, 16-QAM or 4-PSK is 4 or 2
% @iter_num: number of iterations, usually select 10
function x_hat = EP(y,H,noiseLevel,modType,mod_size,iter_num)
damp_mes = 0.9;

% ModulatedSym =qammod(0:2^mod_size-1,2^mod_size, 'UnitAveragePower',true)'; 
% real_mod_Sym = imag(ModulatedSym(1:mod_size));

%  x = [real(x); imag(x)] ;
%                 y = [real(y); imag(y)] ;
%                         H = [real(H) -imag(H); imag(H) real(H) ];   
                        noiseLevel = noiseLevel/2 ; 
						Ht = H' ; 
						var_Noise = noiseLevel;
						[~, U_iter] = size(H);
						gamma = zeros(U_iter,1);
						lam = ones(U_iter,1);
						lamda = 2*lam;
						t= 0;
                        gxIn = ConstellationEstimIn_Real(modType,mod_size, 0, 0.5) ;
%                         gxIn = ConstellationEstimIn(modType,mod_size, 0, 1) ;


while (t < iter_num) 
            
%                              %% Find Zigma and Mu    
                        zigma_A  = inv ( Ht * H + var_Noise * diag(lamda));
                        mu_A = zigma_A *(Ht*y + var_Noise*gamma);
        				
						zigma_A_B = var_Noise * real(diag(zigma_A));
						mu_A_B = mu_A;
%                             %% Find Cavity Distribution
						var_Cav = max(zigma_A_B ./ (1-zigma_A_B .* lamda), eps) ;
						mu_Cav = var_Cav  .* (mu_A_B ./ zigma_A_B - gamma) ;
%                             %% Exp Family Mean and Variance
						[mu_B, zigma_B] = gxIn.estim(mu_Cav, var_Cav) ; 
						zigma_B = max(zigma_B,10^-13);
%                             %% Update Parameters
						lamda_New = (var_Cav -  zigma_B)  ./zigma_B  ./   var_Cav ;
						gamma_New = ( mu_B  .* var_Cav   -  zigma_B  .*  mu_Cav ) ./  zigma_B  ./  var_Cav ;  
% % %                             %% Avoid Negative Index
						negIdx = lamda_New < 0;
						lamda_New(negIdx)  =  lamda(negIdx);
						gamma_New(negIdx)  =  gamma(negIdx);
%                             %% Damping Calculation
						lamda = damping(lamda, lamda_New, damp_mes);
						gamma = damping(gamma, gamma_New, damp_mes);
                        
% cons_real_eq = real(qammod([0:sqrt(2^mod_size):2^mod_size-1]',2^mod_size,'UnitAveragePower',true));                
% x_idx = XTrain_gens(cons_real_eq,x);
% x_hat_idx = XTrain_gens(cons_real_eq,mu_B);
% for idx =1:length(real_mod_Sym)
% prob_s1(:,idx) = normpdf(real_mod_Sym(idx),mu_B,zigma_B);
% end


t = t+1;     

end
% result.error = error;
x_hat = mu_B(1:U_iter/2) + 1j*mu_B(U_iter/2+1:U_iter) ;   
end
