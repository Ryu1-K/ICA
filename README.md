# ICA
SIM.nworker= 15;    % Num. of workers   
SIM.M      = 4;     % Num. of TX antennas
SIM.N      = 16;    % Num. of RX antennas
SIM.Kd     = 512;   % Num. of symbols
SIM.wloop  = 3;     % Num. of Trials
SIM.method = 'ICA';  % Detector (ZF,ICA)
SIM.ml     = 2;     % Modulation level (2:QPSK, 4:16QAM, 6:64QAM)
SIM.EsN0   = [0:1:15]; % Es/N0
