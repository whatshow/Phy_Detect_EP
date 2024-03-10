close all;
clear;
clc;

%% file system
path_folder = "./_tmp/Samples/Tests/";
path_file = path_folder + "test_case_02.mat";
% try to load
if ~exist(path_file, 'file')
    error("The data is not created.")
else
    load(path_file);
end

%% Param Config - Model
SNR_range = 8:2:16;                                        % SNR range
No_range = 10.^(-SNR_range/10);
M = 16;                                                     % M-ary QAM         
tx_num = 12;                                                % Tx antenna number
rx_num = 12;                                                % Rx antenna number

%% Param Config - EP algorithm 
ep_iter_times = 10;                                         % The times to estimate the distribution of  
ep_beta = 0.9;                                              % Convex combination, usually is 0.2
minvar = 1e-13;
early_stop_diff = 1e-4;

%% Simulation
SERs_ep = zeros(1, length(SNR_range));
SERs_ep_es = zeros(1, length(SNR_range));
SERs_alva = zeros(1, length(SNR_range));

syms_all = zeros(tx_num, 1, length(SNR_range), nframe);
syms_es_all = zeros(tx_num, 1, length(SNR_range), nframe);
syms_alva_all = zeros(tx_num, 1, length(SNR_range), nframe);

for idx = 1:length(SNR_range)
    % Get current SNR
    SNR = SNR_range(idx);
    No = No_range(idx);
    fprintf("SNR=%d\n", SNR);
    
    % Prepare the space to store all BERs during 'BER_max_try_times' times
    SERs_ep_tmp = zeros(1, nframe);
    SERs_ep_es_tmp = zeros(1, nframe);
    SERs_alva_tmp = zeros(1, nframe);
    % Try several times to do average on all BERs to avoid fluctuation
    for try_times = 1:nframe
        % to real
        y_real = y_real_all(:, :, idx, try_times);
        H_real = H_real_all(:, :, idx, try_times);
        x= x_all(:, :, idx, try_times);
        
        % EP
        % EP 
        ep = EP(sympool_real, "beta", (1-ep_beta), "epsilon", minvar, "early_stop_min_diff", early_stop_diff);
        syms = ep.detect(y_real,H_real, No/2, "sym_map", true);
        syms = syms(1:length(syms)/2) + 1j*syms(length(syms)/2+1:end);
        % EP - early stop
        ep_es = EP(sympool_real, "beta", (1-ep_beta), "epsilon", minvar, "early_stop", true, "early_stop_min_diff", early_stop_diff);
        syms_es = ep_es.detect(y_real,H_real, No/2, "sym_map", true);
        syms_es = syms_es(1:length(syms_es)/2) + 1j*syms_es(length(syms_es)/2+1:end);
        % EP - alva
        ep_alva_symmap = EP(sympool);
        [syms_alva] = EP_Alva(y_real,H_real, No, 'QAM', log2(M), ep_iter_times);
        syms_alva = ep_alva_symmap.symmap(syms_alva);

        % save detection results
        syms_all(:, :, idx, try_times) = syms;
        syms_es_all(:, :, idx, try_times) = syms_es;
        syms_alva_all(:, :, idx, try_times) = syms_alva;

        % SER cal
        SERs_ep_tmp(try_times) = sum(x - syms > eps)/tx_num;
        SERs_ep_es_tmp(try_times) = sum(x - syms_es > eps)/tx_num;
        SERs_alva_tmp(try_times) = sum(x - syms_alva > eps)/tx_num;
    end
    % mean
    SERs_ep(idx) = mean(SERs_ep_tmp);
    SERs_ep_es(idx) = mean(SERs_ep_es_tmp);
    SERs_alva(idx) = mean(SERs_alva_tmp);
end

%% Save
save(path_file, "SERs_alva", "SERs_ep", "SERs_ep_es", "-append");
save(path_file, "sym_bitnum", "sympool", "sympool_real", "-append");
save(path_file, "syms_all", "syms_es_all", "syms_alva_all", "-append");

%% plot
semilogy(SNR_range, SERs_alva, "-s", "Color", "#D95319", "LineWidth", 4);
hold on;
semilogy(SNR_range, SERs_ep, "--ob", "LineWidth", 2);
hold on;
semilogy(SNR_range, SERs_ep_es, ":>", "Color", "#77AC30", "LineWidth", 2);
hold off;
grid on;
xlabel("SNR(dB)");
ylabel("SER");
ylim([min(SERs_alva), 1]);
xlim([min(SNR_range), max(SNR_range)]);
legend('EP-Alva', 'EP', 'EP-EarlyStop');
