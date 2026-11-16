function report = validate_simulation(out)
%VALIDATE_SIMULATION 默认已标定配置的验收；关闭标定会使相干验收预期失败。
% 核验物理符号、噪声、理想TDM补偿、BP相干性和默认场景分辨率。
% validate_simulation(out) 复用主仿真结果；无参数时运行主仿真后验证。
if nargin<1, out=main_tdm_mimo_sar(struct('figureVisible','off')); end
cfg=out.cfg;
% 1. 直接 TX*conj(RX) 与解析 IF 一致。只以双精度数值误差作为容差。
t=(0:63)'/cfg.fs; tau=2*5.031/cfg.c;
tx=exp(1i*2*pi*(cfg.fc*t+0.5*cfg.slope*t.^2));
rx=exp(1i*2*pi*(cfg.fc*(t-tau)+0.5*cfg.slope*(t-tau).^2));
analytic=exp(1i*2*pi*(cfg.slope*tau*t+cfg.fc*tau-0.5*cfg.slope*tau^2));
dechirpError=max(abs(tx.*conj(rx)-analytic));
assert(dechirpError<1e-8,'去斜解析模型符号错误。');
wave=fmcw_waveform_demo(cfg,out.geom);
assert(wave.maxMixingError<1e-8,'显式发射/接收复包络混频结果不一致。');
% 2. FFT中心化 + RVP补偿：非整bin点的相干恢复。
w=0.5-0.5*cos(2*pi*(0:63)'/63); nfft=1024;
freq=(0:nfft-1)'*cfg.fs/nfft; tc=63/(2*cfg.fs);
yc=fft(analytic.*w,nfft)/sum(w).*exp(1i*2*pi*freq*tc);
value=interp1(freq,yc,cfg.slope*tau,'linear')*exp(-1i*2*pi*(cfg.fc+cfg.slope*tc)*tau+1i*pi*cfg.slope*tau^2);
assert(abs(angle(value))<1e-6 && abs(value)>0.99,'距离FFT插值或BP匹配相位错误。');
% 3. 单点几何BP：不依赖多目标成像肉眼判断的真实相干验收。
single=cfg; single.nLoops=17; single.nFrames=3; single.snrDb=Inf;
single.velocity=[1,1]; % 加入明确的径向运动，使无补偿对照具有可测相位变化
single.targets=[0.03,5.025]; single.reflectivity=1;
[a,g]=simulate_echo(single); [r,ra]=range_compress(a,single);
grid=struct('X',single.targets(1),'Y',single.targets(2),'mask',true);
[im,st]=bp_image(r,ra,g,grid,single,'correct');
[imBad,~]=bp_image(r,ra,g,grid,single,'no_phase');
assert(st.coherence>0.999,'无噪声单目标BP相干系数不足。');
assert(abs(im)>1.5*abs(imBad),'单点补偿增益不足。');
% 4. AWGN功率必须满足设定，且TDM补偿后单目标角度回到真值。
if isfinite(cfg.snrDb), assert(abs(out.echoInfo.measuredSnrDb-cfg.snrDb)<0.1,'实测SNR偏差过大。'); end
demo=out.tdmDemo;
assert(abs(demo.estimatedAfter-demo.trueAngle)<0.1,'TDM补偿后的角度不正确。');
assert(abs(demo.estimatedAfter-demo.trueAngle)<abs(demo.estimatedBefore-demo.trueAngle),'TDM补偿未改善测角。');
% 5. 检查真值栅格在观测生成的ROI内。此处才读取目标真值作为验收。
covered=false(size(cfg.targets,1),1);
for k=1:numel(covered)
    [~,ix]=min(abs(out.mimo.x-cfg.targets(k,1))); [~,iy]=min(abs(out.mimo.y-cfg.targets(k,2)));
    covered(k)=out.mimo.mask(iy,ix);
end
assert(all(covered),'ROI漏掉目标：需检查噪声水平或扩展ROI；不能用真值修补检测。');
% 近邻分辨验收仅对默认足够长孔径的场景执行。
peakPositions=[]; valleyDb=NaN;
defaultNearPoints=[-5*sind(0.5),5*cosd(0.5);5*sind(0.5),5*cosd(0.5)];
if size(cfg.targets,1)>=2 && out.echoInfo.apertureLength>0.38 && ...
        max(abs(cfg.targets(1:2,:)-defaultNearPoints),[],'all')<1e-8 && ...
        min(out.mimo.x)<=-0.09 && max(out.mimo.x)>=0.09
    [~,iy]=min(abs(out.mimo.y-cfg.targets(1,2)));
    x=out.mimo.x; p=abs(out.image(iy,:));
    left=find(x>=-0.09 & x<=-0.01); right=find(x>=0.01 & x<=0.09);
    [pl,jl]=max(p(left)); [pr,jr]=max(p(right)); il=left(jl); ir=right(jr);
    peakPositions=x([il,ir]);
    assert(max(abs(peakPositions-sort(cfg.targets(1:2,1)).'))<=2*cfg.dx,'近邻峰值定位误差超过两像素。');
    valleyDb=20*log10(min(p(il:ir))/min(pl,pr));
    assert(valleyDb<-3,'两个近邻点之间未形成3dB谷值。');
end
report=struct('passed',true,'dechirpMaxError',dechirpError, ...
    'fftMatchedPhaseErrorRad',angle(value),'fftMatchedAmplitude',abs(value), ...
    'singlePointCoherence',st.coherence,'singlePointGainDb',20*log10(abs(im)/abs(imBad)), ...
    'roiCoversTargets',covered,'nearbyPeakPositionsM',peakPositions,'valleyDb',valleyDb);
disp(report);
save(fullfile(cfg.outputDir,'validation_report.mat'),'report');
fid=fopen(fullfile(cfg.outputDir,'验证结果.txt'),'w','n','UTF-8'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'全部验证通过。\n去斜最大误差 %.3g\nFFT匹配相位误差 %.3g rad\nFFT匹配幅度 %.6f\n单点相干系数 %.9f\n单点相位补偿增益 %.3f dB\n',dechirpError,angle(value),abs(value),st.coherence,report.singlePointGainDb);
fprintf(fid,'ROI覆盖全部目标：%d\n近邻峰位置：',all(covered)); fprintf(fid,' %.4f m',peakPositions);
fprintf(fid,'\n近邻目标间谷值：%.3f dB\n',valleyDb);
end
