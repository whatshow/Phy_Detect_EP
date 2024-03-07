close all;
clear;
clc;
%% file
path_folder = "./_dist/Samples/Tests/";
path_file = path_folder + "test_case_01.mat";
% create the folder if not exist
if ~exist(path_folder, 'dir')
    mkdir(path_folder);
end
% delete the file if exist
if exist(path_file, 'file')
    delete(path_file)
end

%% Param Config - Model
SNR = 16;                                                   % SNR range
No = 10^(-SNR/10);
M = 16;                                                     % M-ary QAM
sym_bitnum = log2(M);                                       % Bit number in 1 M-ary modulation symbol
sympool = qammod([0: M - 1], M, "UnitAveragePower", true);  % The symbol pool to store all possible M-ary modulation symbols
sympool_real = unique(real(sympool));              
tx_num = 12;                                                % Tx antenna number
rx_num = 12;                                                % Rx antenna number
% Param Config - EP algorithm 
ep_iter_times = 10;                                         % The times to estimate the distribution of  
ep_beta = 0.9;                                              % Convex combination, usually is 0.2s
minvar = 1e-13;
early_stop_diff = 1e-4;

%% Gen data
nbits_len = tx_num*sym_bitnum;
nbits = randi([0 1], nbits_len, 1);
% Create symbols
x = qammod(nbits, M,'InputType','bit','UnitAveragePower',true);
H = (randn(tx_num, rx_num) + 1j*randn(tx_num, rx_num))/sqrt(2*tx_num);
noise = sqrt(No/2) * (randn(rx_num,1) + 1j*randn(rx_num,1)) ;
% Through AWGN channel to get y 
y = H*x + noise;
x_real = [real(x);imag(x)];
y_real = [real(y);imag(y)];
H_real = [real(H), -imag(H); imag(H), real(H)];


%% save
save(path_file, "M", "sympool", "sympool_real");
save(path_file, "y", "H", "x", "nbits", "No", "y_real", "H_real", "x_real", "-append");
save(path_file, "ep_iter_times", "ep_beta", "minvar", "early_stop_diff", "-append");