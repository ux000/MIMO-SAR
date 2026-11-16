function wave = fmcw_waveform_demo(cfg,geom)
%FMCW_WAVEFORM_DEMO 显式构造发射、延迟接收、混频去斜三种波形。
% 射频形式：sTX=exp(j*2*pi*(fc*t+S*t^2/2))，sRX=sTX(t-tau)。
% 本函数用复包络表示，保留接收载频延迟 exp(-j*2*pi*fc*tau)。
% 展示 chirp 复包络需高采样率；不能用 4MHz 的 IF ADC 直接采样336MHz chirp。
rt=hypot(cfg.targets(1,1)-geom.txOffset(1),cfg.targets(1,2));
rr=hypot(cfg.targets(1,1)-geom.rxOffset(1),cfg.targets(1,2));
tau=(rt+rr)/cfg.c;
duration=cfg.nSamples/cfg.fs;
fsWave=4*cfg.slope*duration; % 大于复包络最高频率的两倍，独立于 IF ADC
t=(0:floor(duration*fsWave)-1)'/fsWave;
tx=exp(1i*pi*cfg.slope*t.^2);
rx=exp(-1i*2*pi*cfg.fc*tau+1i*pi*cfg.slope*(t-tau).^2);
rx(t<tau)=0;              % 第一束回波在传播延迟后到达
mixed=tx.*conj(rx);
analytic=exp(1i*2*pi*(cfg.slope*tau*t+cfg.fc*tau-0.5*cfg.slope*tau^2));
analytic(t<tau)=0;
wave=struct('time',t,'txBaseband',tx,'rxBaseband',rx,'ifSignal',mixed, ...
    'tau',tau,'sampleRate',fsWave,'beatHz',cfg.slope*tau, ...
    'maxMixingError',max(abs(mixed-analytic)));
end
