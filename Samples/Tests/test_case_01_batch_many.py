import numpy as np
import sys
import scipy.io
import os
sys.path.append("..");
from EP import EP
from Samples.Alva.Utils.EP_Alva import EP_Alva;

# batch_size
batch_size = 10;

# file
project_name = "phy_detect_ep";
path_folder = os.path.abspath(os.path.dirname(__file__)).lower();
path_folder = path_folder[:path_folder.find(project_name)+len(project_name)];
path_file = os.path.normpath(path_folder+"/_dist/Samples/Tests/test_case_01.mat");

# load matlab data
try:
    matlab_data = scipy.io.loadmat(path_file);
except FileNotFoundError:
    raise Exception("You have to run `Samples/Tests/test_case_01_gen_data.m` to create data.");
sympool_real = matlab_data["sympool_real"];
beta = 1 - matlab_data["ep_beta"];
epsilon = matlab_data["minvar"].squeeze();
early_stop_min_diff = matlab_data["early_stop_diff"];
y_real = np.tile(matlab_data["y_real"].squeeze(), (batch_size, 1));
H_real = np.tile(matlab_data["H_real"], (batch_size, 1, 1));
No = matlab_data["No"].squeeze();
try:
    syms_alva = np.tile(matlab_data["syms_alva"].squeeze(), (batch_size, 1));
    syms_matlab = np.tile(matlab_data["symbols"].squeeze(), (batch_size, 1));
except KeyError:
    raise Exception("You have to run `Samples/Tests/test_case_01.m` to retrieve estimated symbols in Matlab.");
    
# EP detection
ep = EP(sympool_real, beta=beta, epsilon=epsilon, early_stop_min_diff=early_stop_min_diff, batch_size=batch_size);
syms_ep = ep.detect(y_real, H_real, No/2);
syms_ep = syms_ep[:, :int(syms_ep.shape[-1]/2)] + 1j*syms_ep[:, int(syms_ep.shape[-1]/2):];
# EP detection - Alva's python code
syms_ep_alva = EP_Alva(sympool_real.squeeze(), y_real.shape[-1], y_real, H_real, No/2, 10);
syms_ep_alva = syms_ep_alva[:, :int(syms_ep_alva.shape[-1]/2)] + 1j*syms_ep_alva[:, int(syms_ep_alva.shape[-1]/2):];

# difference
syms_diff0 = abs(syms_ep_alva - syms_alva);
print("The difference between Alva's EP(python) and Alva's EP(matlab) is %e"%(np.sum(syms_diff0, axis=None)/batch_size));
syms_diff1 = abs(syms_ep - syms_alva);
print("The difference between EP(python) and Alva's EP(matlab) is %e"%(np.sum(syms_diff1, axis=None)/batch_size));
syms_diff2 = abs(syms_ep - syms_matlab);
print("The difference between EP(python) and EP(matlab) is %e"%(np.sum(syms_diff2, axis=None)/batch_size));
syms_diff3 = abs(syms_ep - syms_ep_alva);
print("The difference between EP(python) and Alva's EP(python) is %e"%(np.sum(syms_diff3, axis=None)/batch_size));
