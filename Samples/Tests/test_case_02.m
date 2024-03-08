close all;
clear;
clc;

%% file system
path_folder = "./_tmp/Samples/Tests/";
path_file = path_folder + "test_case_02.mat";

%% Param Config - Model
SNR_range = 8:2:16;                                        % SNR range
No_range = 10.^(-SNR_range/10);
M = 16;                                                     % M-ary QAM
sym_bitnum = log2(M);                                       % Bit number in 1 M-ary modulation symbol
sympool = qammod([0: M - 1], M, "UnitAveragePower", true);  % The symbol pool to store all possible M-ary modulation symbols
sympool_real = 1/sqrt((2/3)*(M-1))*pammod([0:(2^(log2(M)/2))-1],(2^(log2(M)/2)));              
tx_num = 12;                                                % Tx antenna number
rx_num = 12;                                                % Rx antenna number
nframe = 1e5;                                               % The times to calculate BER to get the mean BER

%% Param Config - EP algorithm 
ep_iter_times = 10;                                         % The times to estimate the distribution of  
ep_beta = 0.9;                                              % Convex combination, usually is 0.2
minvar = 1e-13;
early_stop_diff = 1e-4;

%% Simulation
SERs_ep = zeros(1, length(SNR_range));
SERs_ep_es = zeros(1, length(SNR_range));
SERs_alva = zeros(1, length(SNR_range));

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
    parfor try_times = 1:nframe
        % nbits
        nbits_len = tx_num*sym_bitnum;
        nbits = randi([0 1], nbits_len, 1);
        % Create symbols
        x = qammod(nbits, M,'InputType','bit','UnitAveragePower',true);
        % Channel 
        %H = 1/sqrt(2*tx_num)*randn(tx_num, rx_num) + 1/sqrt(2)*randn(tx_num, rx_num)*1j;
        H = (randn(tx_num, rx_num) + 1j*randn(tx_num, rx_num))/sqrt(2*tx_num) ;
        % Noise Creation
        noiseLevel = 10^(-SNR/10);
        noise = sqrt(noiseLevel/2) * (randn(rx_num,1) + 1j*randn(rx_num,1)) ;
        %noise = (randn(rx_num, 1) + randn(rx_num, 1)*1j)*sqrt(noiseLv/2);
        % Through AWGN channel to get y 
        y = H*x + noise;
        % to real
        y_real = [real(y);imag(y)];
        H_real = [real(H), -imag(H); imag(H), real(H)];
        
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
        [syms_alva] = EP_Alva(y_real,H_real,noiseLevel, 'QAM', log2(M), ep_iter_times);
        syms_alva = ep_alva_symmap.symmap(syms_alva);

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
save(path_file, "SERs_alva", "SERs_ep", "SERs_ep_es");

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
