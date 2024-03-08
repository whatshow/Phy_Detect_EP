import numpy as np
import sys
import scipy.io
import os
sys.path.append("..");
from EP import EP

# file
project_name = "phy_detect_ep";
path_folder = os.path.abspath(os.path.dirname(__file__)).lower();
path_folder = path_folder[:path_folder.find(project_name)+len(project_name)];
path_file = os.path.normpath(path_folder+"/_tmp/Samples/Tests/test_case_02.mat");

# load matlab data
try:
    matlab_data = scipy.io.loadmat(path_file);
except FileNotFoundError:
    raise Exception("You have to run `Samples/Tests/test_case_01_gen_data.m` to create data.");
SERs_alva = matlab_data["SERs_alva"];
SERs_ep = matlab_data["SERs_ep"];
SERs_ep_es = matlab_data["SERs_ep_es"];