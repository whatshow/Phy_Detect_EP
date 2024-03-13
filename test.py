import os
import shutil

# configuration
path_cur = os.getcwd();
path_dist = "_dist/";
path_dist_pkg = path_dist + "whatshow_phy_detect_ep/";
file_ep = "EP.py"
file_init = "build-init.py";
file_init_dsc = "__init__.py";
file_setup = "build-setup.py";
file_setup_dsc = "setup.py";

# build & upload
os.chdir(path_dist);
os.system("twine upload dist/*");
os.chdir(path_cur);
