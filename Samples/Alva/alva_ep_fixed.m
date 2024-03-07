close all;
clear;
clc;

%% file
path_folder = "./_dist/Samples/Alva/";
path_file = path_folder + "test_ep_fixed.mat";
% create the folder if not exist
if ~exist(path_file, 'file')
    error("The data is not created.")
else
    load(path_file);
end

%% detect
[syms] = Detection_EP(sympool_real, H_real, y_real, No/2, ep_iter_times, "Beta", ep_beta, "MinVariance", 1e-13);
syms = [syms(1:length(syms)/2) + 1j*syms(length(syms)/2+1:end)];
[syms_alva] = EP_Alva(y_real,H_real, No, 'QAM', log2(M), ep_iter_times);
% check difference
syms_diff = abs(syms - syms_alva);
fprintf("The difference is %.16f\n", sum(syms_diff));