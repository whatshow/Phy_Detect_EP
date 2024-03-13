close all;
clear;
clc;

%% file system
path_folder = "./_tmp/Samples/Tests/";
path_file = path_folder + "test_case_02.mat";
% create the folder if not exist
if ~exist(path_folder, 'dir')
    mkdir(path_folder);
end
% delete the file if exist
if exist(path_file, 'file')
    delete(path_file)
end

%% Param Config - Model
SNR_range = 8:2:16;                                        % SNR range
No_range = 10.^(-SNR_range/10);
M = 16;                                                     % M-ary QAM
sym_bitnum = log2(M);                                       % Bit number in 1 M-ary modulation symbol
sympool = qammod([0: M - 1], M, "UnitAveragePower", true);  % The symbol pool to store all possible M-ary modulation symbols
sympool_real = 1/sqrt((2/3)*(M-1))*pammod([0:(2^(log2(M)/2))-1],(2^(log2(M)/2)));              
tx_num = 12;                                                % Tx antenna number
rx_num = 12;                                                % Rx antenna number
nframe = 100;                                               % The times to calculate BER to get the mean BER

%% Param Config - EP algorithm 
ep_iter_times = 10;                                         % The times to estimate the distribution of  
ep_beta = 0.9;                                              % Convex combination, usually is 0.2
minvar = 1e-13;
early_stop_diff = 1e-4;

%% Simulation
SERs_ep = zeros(1, length(SNR_range));
SERs_ep_es = zeros(1, length(SNR_range));
SERs_alva = zeros(1, length(SNR_range));

x_all = zeros(tx_num, 1, nframe, length(SNR_range));
H_real_all = zeros(rx_num*2, tx_num*2, nframe, length(SNR_range));
y_real_all = zeros(tx_num*2, 1, nframe, length(SNR_range));

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
        % save
        x_all(:, :, try_times, idx) = x;
        H_real_all(:, :, try_times, idx) = H_real;
        y_real_all(:, :, try_times, idx) = y_real;
    end
end

%% Save
save(path_file, "x_all", "H_real_all", "y_real_all");
save(path_file, "nframe", "sym_bitnum", "sympool", "sympool_real", "-append");
