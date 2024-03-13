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
file_readme = "README.md";

# clear the distribution
if os.path.exists(path_dist):
    shutil.rmtree(path_dist);

# create the distribution folder
os.makedirs(path_dist_pkg);    
    
# copy files to the distribution
shutil.copyfile(file_ep, path_dist_pkg + file_ep);
shutil.copyfile(file_init, path_dist_pkg + file_init_dsc);
shutil.copyfile(file_setup, path_dist + file_setup_dsc);
shutil.copyfile(file_readme, path_dist + file_readme);


# build & upload
os.chdir(path_dist);
os.system("python " + file_setup_dsc + " sdist bdist_wheel");
os.system("twine upload dist/*");
os.chdir(path_cur);
