function RES = main_task(SIM,idx_worker)

%% 乱数のシードを設定（worker毎に異なる値）
s = RandStream('mt19937ar','Seed', 1234*idx_worker);
RandStream.setGlobalStream(s);

%% パイロット系列生成
switch(SIM.method)
    case 'ZF' % 直交パイロット
        UW = hadamard(SIM.M);
        UW = UW(1:SIM.M,:);
        TX.UW = UW(1:SIM.M,:);
    case 'ICA' % 非直交パイロット
        UW = 2*de2bi(0:SIM.M-1,'left-msb')-1;
        TX.UW = [ones(SIM.M,1) UW];
end
SIM.Kp = size(TX.UW,2);

%% シミュレーション
for idx_En = 1:length(SIM.EsN0)
    CH.N0 = 10^(-SIM.EsN0(idx_En)/10)*SIM.M; %雑音分散
    CH.sigma = sqrt(CH.N0/2);
    err.noe = zeros(3,1); err.nos = zeros(3,1);
    err.sqerr = 0;
    for idx_loop = 1:ceil(SIM.nloop/SIM.nworker)
        %%% TX %%%
        % ビット生成
        TX.alp = randi(SIM.Q,SIM.M,SIM.Kd)-1;
        TX.bit = de2bi(TX.alp(:),SIM.ml,'left-msb');
        % 変調
        TX.sym = qammod(TX.alp,SIM.Q,'UnitAveragePower',true);
        % パイロット付加
        TX.sig = [TX.UW TX.sym];
        
        %%% CHANNEL %%%
        % フェーディング係数生成
        CH.H = (randn(SIM.N,SIM.M)+randn(SIM.N,SIM.M)*1i)/sqrt(2);
        % 雑音生成
        CH.n = CH.sigma*(randn(SIM.N,SIM.Kd+SIM.Kp)+randn(SIM.N,SIM.Kd+SIM.Kp)*1i);
        % 一様フェージング通信路
        RX.sig = CH.H * TX.sig + CH.n;
        
        %%% RX %%%
        % 検出器
        switch(SIM.method)
            case 'ZF'
                CH.H_ = RX.sig(:,1:SIM.Kp)*TX.UW'/SIM.Kp;
                y = pinv(CH.H_)*RX.sig(:,SIM.Kp+1:end);
            case 'ICA'
                [RX.sym_, ~] = fICA(RX.sig, CH.N0, SIM.M);
                
                %　順序整理
                C = real( (RX.sym_(:,2:SIM.Kp)./ RX.sym_(:,1)) * TX.UW(:,2:end)');
                [~,index] = max(C);
                T = zeros(SIM.M);
                T(index + [0:SIM.M:(SIM.M-1)*SIM.M]) = 1;
                RX.sym_ = T'*RX.sym_;

                % 位相回転補償
                H_ = mean(RX.sym_(:,1:SIM.Kp).*TX.UW,2);
                y= RX.sym_(:,SIM.Kp+1:end)./H_;
        end
        
        % 復調器
        RX.alp = qamdemod(y,2^SIM.ml,'UnitAveragePower',true);
        RX.bit = de2bi(RX.alp(:),SIM.ml,'left-msb');
        
        % 通信路推定
        RX.UW = qammod(RX.alp,SIM.Q,'UnitAveragePower',true);
        CH.H_hat = RX.sig*pinv([TX.UW RX.UW]);
        
        % 最尤検出
        % 自分で作る。
                
        %%% ERR %%%
        noe_ins = sum(sum(RX.bit~=TX.bit));
        err.noe(1) = err.noe(1)+noe_ins;
        err.noe(2) = err.noe(2)+ sum(sum(RX.alp ~= TX.alp));
        err.noe(3) = err.noe(3)+(noe_ins~=0); 
        err.nos(1) = err.nos(1)+SIM.M*SIM.ml*SIM.Kd;
        err.nos(2) = err.nos(2)+SIM.M*SIM.Kd;
        err.nos(3) = err.nos(3)+1;
        err.sqerr  = err.sqerr + mean(mean(abs(CH.H_hat - CH.H).^2));
    end
    RES(idx_En, :) = [err.noe; err.nos; err.sqerr];
end