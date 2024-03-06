import numpy as np

'''
EP detection
<INPUTS>
@y_noise_real:            the received y
@H_real:                  the channel
@sigma2:                  the noise power (in the real equivalent domain)
@num_iter:                the number of iteration for EP
'''
def EP_Alva(constellation, user_num, x_real,y_noise_real,H_real,sigma2,num_iter):
    lamda_init = np.ones((H_real.shape[0], user_num))*2
    gamma_init = np.zeros((H_real.shape[0], user_num));
    sigma2 = np.mean(sigma2)
    H = H_real
    y = y_noise_real
    constellation_expanded = np.expand_dims(constellation, axis=1)
    constellation_expanded= np.repeat(constellation_expanded[None,...],H.shape[0],axis=0)
    

    def calculate_mean_var(pyx, constellation_expanded):
        constellation_expanded_transpose = np.repeat(constellation_expanded.transpose(0,2,1), user_num, axis=1)
        mean = np.matmul(pyx, constellation_expanded)
        var = np.square(np.abs(constellation_expanded_transpose - mean))
        var = np.multiply(pyx, var) 
        var = np.sum(var, axis=2)
        
        return np.squeeze(mean), var
    
    def calculate_pyx( mean, var, constellation_expanded):
        constellation_expanded_transpose = np.repeat(constellation_expanded.transpose(0,2,1), user_num, axis=1)
        arg_1 = np.square(np.abs(constellation_expanded_transpose - np.expand_dims(mean,2)))
        log_pyx = (-1 * arg_1)/(2*np.expand_dims(var,2))
        log_pyx = log_pyx - np.expand_dims(np.max(log_pyx,2),2)
        p_y_x = np.exp(log_pyx)
        p_y_x = p_y_x/(np.expand_dims(np.sum(p_y_x, axis=2),2) + np.finfo(float).eps)
        
        return p_y_x

    def LMMSE( H, y, sigma2, lamda, gamma):
        HtH = np.matmul(np.transpose(H,[0,2,1]), H)
        Hty = np.squeeze(np.matmul(np.transpose(H,[0,2,1]), np.expand_dims(y,2)))
        diag_lamda = np.zeros((HtH.shape[0],user_num,user_num))
        np.einsum('ijj->ij',diag_lamda)[...] = lamda
        var = np.linalg.inv(HtH + diag_lamda * sigma2)
        mean = (Hty) + gamma* sigma2
        mean = np.matmul(var,np.expand_dims(mean,2))
        var = var* sigma2
        return mean, var

    lamda = lamda_init
    gamma = gamma_init

    for iteration in range(num_iter):
                              
        mean_mmse, var_mmse = LMMSE(H, y, sigma2, lamda, gamma)

        # Calculating mean and variance of P_y_x
        
        diag_mmse=np.diagonal(var_mmse, axis1=1, axis2=2)
        var_ab = (diag_mmse/ (1 - diag_mmse*lamda )) +np.finfo(float).eps
        
        # var_ab = 1 / (1/np.diagonal(var_mmse, axis1=1, axis2=2) - lamda )
        mean_ab = (np.squeeze(mean_mmse)/np.diagonal(var_mmse, axis1=1, axis2=2) - gamma) 
        mean_ab = var_ab * mean_ab

        # Calculating P_y_x
        p_y_x_ab = calculate_pyx (mean_ab, var_ab, constellation_expanded)

        # Calculating mean and variance of \hat{P}_x_y
        mean_b, var_b = calculate_mean_var(p_y_x_ab, constellation_expanded)
        var_b = np.clip(var_b, 1e-13, None)

        # Calculating new lamda and gamma
        lamda_new = ((var_ab-var_b) / var_b )/ var_ab
        
        # lamda_new = 1/ var_b - 1/var_ab
        gamma_new = mean_b /var_b - mean_ab/ var_ab

        # Avoiding negative lamda and gamma
        if np.any(lamda_new < 0):
            indices = np.where(lamda_new<0)
            lamda_new[indices]=lamda[indices]
            gamma_new[indices]=gamma[indices]

        # Appliying updating weight
        lamda = lamda*0.9 + lamda_new*0.1
        gamma = gamma*0.9 + gamma_new*0.1
    
    return mean_b;