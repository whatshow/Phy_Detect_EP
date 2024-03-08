close all;
clear;
clc;
%%
path_folder = "./_tmp/Samples/Tests/";
path_file = path_folder + "test_case_01.mat";
% create the folder if not exist
if ~exist(path_file, 'file')
    error("The data is not created.")
else
    load(path_file);
end

%% detect
% detect - By Xinwei
[syms_xin] = Detection_EP(sympool_real, H_real, y_real, No/2, ep_iter_times, "Beta", ep_beta, "MinVariance", minvar);
syms_xin = [syms_xin(1:length(syms_xin)/2) + 1j*syms_xin(length(syms_xin)/2+1:end)];
% detect - by Alva
[syms_alva] = EP_Alva(y_real,H_real, No, 'QAM', log2(M), ep_iter_times);
% detect - EP
ep = EP(sympool_real, "beta", (1-ep_beta), "epsilon", minvar, "early_stop_min_diff", early_stop_diff);
symbols = ep.detect(y_real,H_real, No/2);
symbols = [symbols(1:length(symbols)/2) + 1j*symbols(length(symbols)/2+1:end)];

% check difference
syms_diff = abs(syms_xin - syms_alva);
fprintf("Original Code: the difference is %.16f\n", sum(syms_diff));
fprintf("New Code: (1) the difference is %.16f\n", sum(abs(syms_xin - symbols)));
fprintf("          (2) the difference is %.16f\n", sum(abs(syms_alva - symbols)));

%% save
save(path_file, "syms_alva", "-append");
save(path_file, "symbols", "-append");