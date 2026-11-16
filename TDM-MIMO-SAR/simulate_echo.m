function [adc, geom, info] = simulate_echo(cfg)
%SIMULATE_ECHO 按真实 TDM 时间和双站路程生成去斜复基带回波。
% adc(快时间, 每Tx慢时间循环, Tx, Rx, 帧)，不是把各 Tx 当作同时发射。
% stop-and-go：单 chirp 内冻结位置；不同 chirp 的位置严格按时间更新。
% 目标静止、单次散射，忽略多径、遮挡及 chirp 内运动。
rng(cfg.seed,'twister');
lambda = cfg.c/cfg.fc;
geom.txOffset = ((0:cfg.nTx-1)-(cfg.nTx-1)/2)*cfg.nRx*lambda/2;
geom.rxOffset = ((0:cfg.nRx-1)-(cfg.nRx-1)/2)*lambda/2;
geom.loopTime = cfg.nTx*cfg.slotTime;
assert(cfg.nLoops*geom.loopTime <= cfg.frameTime,'一帧有效 chirp 长度超过帧周期。');
assert(cfg.nSamples/cfg.fs <= cfg.slotTime,'ADC 采样时长超过 Tx 时间槽。');
assert(numel(cfg.channelPhaseDeg)==cfg.nTx*cfg.nRx,'通道相位参数数目不正确。');
assert(size(cfg.targets,1)==numel(cfg.reflectivity),'目标与散射系数数目不一致。');
tFast = (0:cfg.nSamples-1)'/cfg.fs;
geom.time = zeros(cfg.nLoops,cfg.nTx,cfg.nFrames);
adc = complex(zeros(cfg.nSamples,cfg.nLoops,cfg.nTx,cfg.nRx,cfg.nFrames));
for f = 1:cfg.nFrames
    for m = 1:cfg.nTx
        t = (f-1)*cfg.frameTime+(0:cfg.nLoops-1)*geom.loopTime+(m-1)*cfg.slotTime;
        geom.time(:,m,f) = t(:);
        px = cfg.velocity(1)*t; py = cfg.velocity(2)*t;
        for r = 1:cfg.nRx
            z = complex(zeros(cfg.nSamples,cfg.nLoops));
            for k = 1:size(cfg.targets,1)
                rt = hypot(cfg.targets(k,1)-px-geom.txOffset(m),cfg.targets(k,2)-py);
                rr = hypot(cfg.targets(k,1)-px-geom.rxOffset(r),cfg.targets(k,2)-py);
                tau = (rt+rr)/cfg.c; % Tx -> 目标 -> Rx，精确双站距离
                amp = cfg.reflectivity(k);
                if cfg.rangeAttenuation
                    amp = amp*cfg.referenceRange^2./(rt.*rr);
                end
                % TX*conj(RX)：正号去斜。保留 residual video phase (RVP)。
                phase = 2*pi*(cfg.slope*tFast.*tau+cfg.fc*tau-0.5*cfg.slope*tau.^2);
                z = z+amp.*exp(1i*phase);
            end
            ch = (m-1)*cfg.nRx+r;
            adc(:,:,m,r,f) = z*exp(1i*pi/180*cfg.channelPhaseDeg(ch));
        end
    end
end
info.signalPower = mean(abs(adc(:)).^2);
info.noisePower = info.signalPower/10^(cfg.snrDb/10);
if isinf(cfg.snrDb) && cfg.snrDb>0, info.noisePower = 0; end
noise = sqrt(info.noisePower/2)*(randn(size(adc))+1i*randn(size(adc)));
info.measuredNoisePower = mean(abs(noise(:)).^2);
info.measuredSnrDb = 10*log10(info.signalPower/max(info.measuredNoisePower,realmin));
adc = adc+noise;
info.bandwidth = cfg.slope*cfg.nSamples/cfg.fs;
info.rangeResolution = cfg.c/(2*info.bandwidth);
info.apertureLength = norm(cfg.velocity)*(max(geom.time(:))-min(geom.time(:)));
info.maxSingleSlotDisplacement = norm(cfg.velocity)*cfg.slotTime;
end
