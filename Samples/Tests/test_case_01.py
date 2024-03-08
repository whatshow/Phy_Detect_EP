import numpy as np
import sys
import scipy.io
import os
sys.path.append("..");
from EP import EP

project_name = "phy_detect_ep";
path_folder = os.path.abspath(os.path.dirname(__file__)).lower();
path_folder = path_folder[:path_folder.find(project_name)+len(project_name)];
path_file = os.path.normpath(path_folder+"/_tmp/Samples/Tests/test_case_01.mat");

# load matlab data
try:
    matlab_data = scipy.io.loadmat(path_file);
except FileNotFoundError:
    raise Exception("You have to run `Samples/Tests/test_case_01_gen_data.m` to create data.");
sympool_real = matlab_data["sympool_real"];
beta = 1 - matlab_data["ep_beta"];
epsilon = matlab_data["minvar"].squeeze();
early_stop_min_diff = matlab_data["early_stop_diff"];
y_real = matlab_data["y_real"].squeeze();
H_real = matlab_data["H_real"];
No = matlab_data["No"].squeeze();
try:
    syms_alva = matlab_data["syms_alva"].squeeze();
    syms_matlab = matlab_data["symbols"].squeeze();
except KeyError:
    raise Exception("You have to run `Samples/Tests/test_case_01.m` to retrieve estimated symbols in Matlab.");
    
# EP detection
ep = EP(sympool_real, beta=beta, epsilon=epsilon, early_stop_min_diff=early_stop_min_diff);
syms_ep = ep.detect(y_real, H_real, No/2);
syms_ep = syms_ep[:int(len(syms_ep)/2)] + 1j*syms_ep[int(len(syms_ep)/2):];


# difference
syms_diff1 = abs(syms_ep - syms_alva);
print("The difference between EP(python) and Alva's original code is %e"%sum(syms_diff1));
syms_diff2 = abs(syms_ep - syms_matlab);
print("The difference between EP(python) and EP(matlab) is %e"%sum(syms_diff2));