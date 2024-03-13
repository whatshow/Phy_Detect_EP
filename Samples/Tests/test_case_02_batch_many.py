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
path_file_py_batch_0 = os.path.normpath(path_folder+"/_tmp/Samples/Tests/test_case_02_py_batch_0.mat");

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
# load batch_0 data
try:
    matlab_data_batch_0 = scipy.io.loadmat(path_file_py_batch_0);
except FileNotFoundError:
    raise Exception("You have to run `Samples/Tests/test_case_02.py` to create data.");
sr_batch_0_SERs_count_ep = matlab_data_batch_0["sr_batch_0_SERs_count_ep"];
sr_batch_0_SERs_count_ep_es = matlab_data_batch_0["sr_batch_0_SERs_count_ep_es"];

# batch
# batch_size = multiprocessing.cpu_count()*64;
batch_size = 50;

# Param Config - Model
SNR_range = np.arange(8, 17, 2);                            # SNR range
No_range = 10**(-SNR_range/10);
M = 16;                                                     # M-ary QAM         
tx_num = 12;                                                # Tx antenna number
rx_num = 12;                                                # Rx antenna number

# Param Config - EP algorithm 
ep_iter_times = 50;                                         # The times to estimate the distribution of  
ep_beta = 0.9;                                              # Convex combination, usually is 0.2
minvar = 1e-13;
early_stop_diff = 1e-4;

# Simulation
ep = EP(sympool_real, beta=(1-ep_beta), epsilon=minvar, early_stop_min_diff=early_stop_diff, batch_size=batch_size);
ep_es = EP(sympool_real, beta=(1-ep_beta), epsilon=minvar, early_stop=True, early_stop_min_diff=early_stop_diff, batch_size=batch_size);
SERs_ep = np.zeros(len(SNR_range));
SERs_ep_es = np.zeros(len(SNR_range));
SERs_alva = np.zeros(len(SNR_range));
for idx in range(len(SNR_range)):
    # Get current SNR
    SNR = SNR_range[idx];
    No = No_range[idx];
    
    # Prepare the space to store all BERs during 'BER_max_try_times' times
    SERs_ep_tmp = np.zeros(int(nframe/batch_size));
    SERs_ep_es_tmp = np.zeros(int(nframe/batch_size));
    SERs_alva_tmp = np.zeros(int(nframe/batch_size));
    # Try several times to do average on all BERs to avoid fluctuation
    print('\rSNR={:>2d}: 0%        '.format(SNR), end="");
    for try_times in range(int(nframe/batch_size)):
        print('\rSNR={:>2d}: {:.2%}'.format(SNR, try_times*batch_size/nframe), end="");
        # Create symbols
        y_real = y_real_all[:, :, try_times*batch_size:(try_times*batch_size+batch_size), idx];
        H_real = H_real_all[:, :, try_times*batch_size:(try_times*batch_size+batch_size), idx];
        x= x_all[:, :, try_times*batch_size:(try_times*batch_size+batch_size), idx];
        y_real = np.moveaxis(y_real, -1, 0);
        H_real = np.moveaxis(H_real, -1, 0);
        x = np.moveaxis(x, -1, 0);
        # retrieve results
        syms_mat_ep = syms_mat_ep_all[:, :, try_times*batch_size:(try_times*batch_size+batch_size), idx];
        syms_mat_es = syms_mat_es_all[:, :, try_times*batch_size:(try_times*batch_size+batch_size), idx];
        syms_mat_alva = syms_mat_alva_all[:, :, try_times*batch_size:(try_times*batch_size+batch_size), idx];
        syms_mat_ep = np.moveaxis(syms_mat_ep, -1, 0).squeeze(-1);
        syms_mat_es = np.moveaxis(syms_mat_es, -1, 0).squeeze(-1);
        syms_mat_alva = np.moveaxis(syms_mat_alva, -1, 0).squeeze(-1);
        
        # detect
        # EP
        syms = ep.detect(y_real,H_real, No/2, sym_map=True);
        syms = syms[:, :int(syms.shape[-1]/2)] + 1j*syms[:, int(syms.shape[-1]/2):];
        syms_diff_ep = abs(syms_mat_ep - syms);
        if np.sum(syms_diff_ep, axis=None) > eps*batch_size:
            raise Exception("EP is different in %d:%d"%(idx, try_times));
        # EP - early stop;
        syms_es = ep_es.detect(y_real,H_real, No/2, sym_map=True);
        syms_es = syms_es[:, :int(syms_es.shape[-1]/2)] + 1j*syms_es[:, int(syms_es.shape[-1]/2):];
        syms_diff_ep_es = abs(syms_es - syms_mat_es);
        if np.sum(syms_diff_ep_es, axis=None) > eps*batch_size:
            raise Exception("EP(ES) is different in %d:%d"%(idx, try_times));
        # EP - alva
        # EP detection - Alva's python code
        syms_ep_alva = EP_Alva(sympool_real.squeeze(), y_real.squeeze(-1).shape[-1], y_real.squeeze(-1), H_real, No/2, 10);
        syms_ep_alva = ep.symmap(syms_ep_alva);
        syms_ep_alva = syms_ep_alva[:, :int(syms_ep_alva.shape[-1]/2)] + 1j*syms_ep_alva[:, int(syms_ep_alva.shape[-1]/2):];
        syms_diff_ep_alva = abs(syms_ep_alva - syms_mat_alva);
        if np.sum(syms_diff_ep_alva, axis=None) > eps*batch_size:
            raise Exception("EP(Alva) is different in %d:%d"%(idx, try_times));
        
        # diff detect
        SERs_count_ep = (x.squeeze(-1) - syms > eps).astype(float);
        SERs_count_ep_es = (x.squeeze(-1) - syms_es > eps).astype(float);
        SERs_count_ep_alva = (x.squeeze(-1) - syms_ep_alva > eps).astype(float);
        # diff detect compar
        SERs_count_diff_ep = SERs_count_ep - sr_batch_0_SERs_count_ep[idx, try_times*batch_size:(try_times*batch_size+batch_size)];
        if np.sum(SERs_count_diff_ep, axis=None) > 0:
            raise Exception("Error count (EP) is differenct in %d:%d"%(idx, try_times));
        SERs_count_diff_ep_es  = SERs_count_ep_es - sr_batch_0_SERs_count_ep_es[idx, try_times*batch_size:(try_times*batch_size+batch_size)];
        if np.sum(SERs_count_diff_ep_es, axis=None) > 0:
            raise Exception("Error count (EP-Es) is differenct in %d:%d"%(idx, try_times));
        # SER cal
        SERs_ep_tmp[try_times] = np.sum(SERs_count_ep, axis=None)/float(tx_num*batch_size);
        SERs_ep_es_tmp[try_times] = np.sum(SERs_count_ep_es, axis=None)/float(tx_num*batch_size);
        SERs_alva_tmp[try_times] = np.sum(SERs_count_ep_alva, axis=None)/float(tx_num*batch_size);
    print('\rSNR={:>2d}: 100%       '.format(SNR));
    # mean
    SERs_ep[idx] = mean(SERs_ep_tmp);
    SERs_ep_es[idx] = mean(SERs_ep_es_tmp);
    SERs_alva[idx] = mean(SERs_alva_tmp);
    
# calculate diff
SERs_diff_ep = SERs_ep - SERs_ep_mat;
SERs_diff_ep_es = SERs_ep_es - SERs_ep_es_mat;
SERs_diff_alva = SERs_alva - SERs_alva_mat;

print("EP difference 1 is %e"%np.sum(SERs_diff_ep));
print("EP difference 2 is %e"%np.sum(SERs_diff_ep_es));
print("EP difference 3 is %e"%np.sum(SERs_diff_alva));