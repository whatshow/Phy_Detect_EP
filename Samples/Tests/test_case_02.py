import numpy as np
import multiprocessing
from numpy import sqrt, mean
from numpy.random import randint
from numpy.random import randn
from numpy import real, imag
import sys
import scipy.io
import os
sys.path.append("..");
from EP import EP
from Samples.Alva.Utils.EP_Alva import EP_Alva;
eps = np.finfo(float).eps;

# file
project_name = "phy_detect_ep";
path_folder = os.path.abspath(os.path.dirname(__file__)).lower();
path_folder = path_folder[:path_folder.find(project_name)+len(project_name)];
path_file = os.path.normpath(path_folder+"/_tmp/Samples/Tests/test_case_02.mat");
path_file_py =  os.path.normpath(path_folder+"/_tmp/Samples/Tests/test_case_02_py_batch_0.mat");

# load matlab data
try:
    matlab_data = scipy.io.loadmat(path_file);
except FileNotFoundError:
    raise Exception("You have to run `Samples/Tests/test_case_01_gen_data.m` to create data.");
SERs_alva_mat = matlab_data["SERs_alva"];
SERs_ep_mat = matlab_data["SERs_ep"];
SERs_ep_es_mat = matlab_data["SERs_ep_es"];
sympool = matlab_data["sympool"].squeeze();
sympool_real = matlab_data["sympool_real"].squeeze();
sym_bitnum = matlab_data["sym_bitnum"].squeeze();
nframe = matlab_data["nframe"].squeeze();
x_all = matlab_data["x_all"];
H_real_all = matlab_data["H_real_all"];
y_real_all = matlab_data["y_real_all"];
# detection results
syms_mat_ep_all = matlab_data["syms_all"];
syms_mat_es_all = matlab_data["syms_es_all"];
syms_mat_alva_all = matlab_data["syms_alva_all"];


# Param Config - Model
SNR_range = np.arange(8, 17, 2);                            # SNR range
No_range = 10**(-SNR_range/10);
M = 16;                                                     # M-ary QAM         
tx_num = 12;                                                # Tx antenna number
rx_num = 12;                                                # Rx antenna number

# Param Config - EP algorithm 
ep_iter_times = 10;                                         # The times to estimate the distribution of  
ep_beta = 0.9;                                              # Convex combination, usually is 0.2
minvar = 1e-13;
early_stop_diff = 1e-4;

# Simulation
# Sim result
sr_batch_0_SERs_count_ep = np.zeros((len(SNR_range), nframe, tx_num)); 
sr_batch_0_SERs_count_ep_es = np.zeros((len(SNR_range), nframe, tx_num));
# Sim - parameters
ep = EP(sympool_real, beta=(1-ep_beta), epsilon=minvar, early_stop_min_diff=early_stop_diff);
ep_es = EP(sympool_real, beta=(1-ep_beta), epsilon=minvar, early_stop=True, early_stop_min_diff=early_stop_diff);
SERs_ep = np.zeros(len(SNR_range));
SERs_ep_es = np.zeros(len(SNR_range));
SERs_alva = np.zeros(len(SNR_range));
for idx in range(len(SNR_range)):
    # Get current SNR
    SNR = SNR_range[idx];
    No = No_range[idx];
    
    # Prepare the space to store all BERs during 'BER_max_try_times' times
    SERs_ep_tmp = np.zeros(nframe);
    SERs_ep_es_tmp = np.zeros(nframe);
    SERs_alva_tmp =np.zeros(nframe);
    # Try several times to do average on all BERs to avoid fluctuation
    print('\rSNR={:>2d}: 0%        '.format(SNR), end="");
    for try_times in range(nframe):
        print('\rSNR={:>2d}: {:.2%}'.format(SNR, try_times/nframe), end="");
        # # Create symbols
        # x_idx = randint(M, size=(batch_size, tx_num, 1));
        # x = np.take(sympool, x_idx)
        # # Channel 
        # H = (randn(batch_size, tx_num, rx_num) + 1j*randn(batch_size, tx_num, rx_num))/sqrt(2*tx_num);
        # # Noise Creation
        # noise = sqrt(No/2) * (randn(batch_size, rx_num, 1) + 1j*randn(batch_size, rx_num, 1));
        # # Through AWGN channel to get y 
        # y = H @ x;
        # # to real
        # x_real = np.concatenate([real(x), imag(x)], axis=-2);
        # y_real = np.concatenate([real(y), imag(y)], axis=-2);
        # H_real = np.concatenate([np.concatenate([real(H), -imag(H)], axis=-1), np.concatenate([imag(H), real(H)], axis=-1)], axis=-2);
        y_real = y_real_all[:, :, try_times, idx];
        H_real = H_real_all[:, :, try_times, idx];
        x= x_all[:, :, try_times, idx];
        # retrieve results
        syms_mat_ep = syms_mat_ep_all[:, :, try_times, idx];
        syms_mat_es = syms_mat_es_all[:, :, try_times, idx];
        syms_mat_alva = syms_mat_alva_all[:, :, try_times, idx];
        syms_mat_ep = syms_mat_ep.squeeze(-1);
        syms_mat_es = syms_mat_es.squeeze(-1);
        syms_mat_alva = syms_mat_alva.squeeze(-1);
        
        # detect
        # EP
        syms = ep.detect(y_real,H_real, No/2, sym_map=True);
        syms = syms[:int(syms.shape[-1]/2)] + 1j*syms[int(syms.shape[-1]/2):];
        syms_diff_ep = abs(syms_mat_ep - syms);
        if np.sum(syms_diff_ep, axis=None) > eps:
            raise Exception("EP is different in %d:%d"%(idx, try_times));
        # EP - early stop
        syms_es = ep_es.detect(y_real,H_real, No/2, sym_map=True);
        syms_es = syms_es[:int(syms_es.shape[-1]/2)] + 1j*syms_es[int(syms_es.shape[-1]/2):];
        syms_diff_ep_es = abs(syms_es - syms_mat_es);
        if np.sum(syms_diff_ep_es, axis=None) > eps:
            raise Exception("EP(ES) is different in %d:%d"%(idx, try_times));

        # diff detect
        SERs_count_ep = (x.squeeze(-1) - syms > eps).astype(int);
        SERs_count_ep_es = (x.squeeze(-1) - syms_es > eps).astype(int);
        # SER cal
        SERs_ep_tmp[try_times] = np.sum(SERs_count_ep, axis=None)/tx_num;
        SERs_ep_es_tmp[try_times] = np.sum(SERs_count_ep_es, axis=None)/tx_num;
        # record
        sr_batch_0_SERs_count_ep[idx, try_times, :] = SERs_count_ep; 
        sr_batch_0_SERs_count_ep_es[idx, try_times, :] = SERs_count_ep_es;
    print('\rSNR={:>2d}: 100%       '.format(SNR));
    # mean
    SERs_ep[idx] = mean(SERs_ep_tmp);
    SERs_ep_es[idx] = mean(SERs_ep_es_tmp);
    
# calculate diff
SERs_ep_diff = SERs_ep - SERs_ep_mat;
SERs_ep_es_diff = SERs_ep_es - SERs_ep_es_mat;

print("EP difference 1 is %e"%np.sum(SERs_ep_diff));
print("EP difference 2 is %e"%np.sum(SERs_ep_es_diff));

# save
scipy.io.savemat(path_file_py, {"sr_batch_0_SERs_count_ep": sr_batch_0_SERs_count_ep, "sr_batch_0_SERs_count_ep_es": sr_batch_0_SERs_count_ep_es});