clear
clc

SIM.nworker= 15;    % Num. of workers PCのスレッド数が上限 (B0001-B0005: 最大8, B0006-B0010: 最大16)     
SIM.M      = 4;     % Num. of TX antennas
SIM.N      = 16;    % Num. of RX antennas
SIM.Kd     = 512;   % Num. of symbols (16QAMの場合には1024以上はほしい)
SIM.wloop  = 3;     % Num. of Trials
SIM.method = 'ICA';  % Detector (ZF,ICA)
SIM.ml     = 2;     % Modulation level (2:QPSK, 4:16QAM, 6:64QAM)
SIM.EsN0   = [0:1:15]; % Es/N0

SIM.nloop  = 10^SIM.wloop;
SIM.Q = 2^SIM.ml;
RES = zeros(length(SIM.EsN0),7);

tic;
if(SIM.nworker==1)
    RES = main_task(SIM,1); %For bug fix
else
    parfor idx_worker = 1:SIM.nworker
        RES_ = main_task(SIM,idx_worker);
        RES = RES +  RES_;
    end
end
toc;

SIM.BER = RES(:,1)./RES(:,4);
SIM.SER = RES(:,2)./RES(:,5);
SIM.FER = RES(:,3)./RES(:,6);
SIM.MSE = RES(:,7)/SIM.nloop;

plot_ber;
fn =  [SIM.method '_' int2str(SIM.M) '_' int2str(SIM.N) '_' int2str(SIM.Kd) '_' int2str(SIM.ml) '_' int2str(SIM.wloop) '.mat'];
save(['DATA\' fn],'SIM')
