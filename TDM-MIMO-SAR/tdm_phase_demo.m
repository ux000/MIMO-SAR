function demo = tdm_phase_demo(cfg,geom)
%TDM_PHASE_DEMO 单目标受控实验：单独放大异步测角误差，避免默认侧视
% 场景小径向速度使 TDM 补偿差异难以观察。该实验不参与主图或ROI。
theta=20; vr=8; fRef=cfg.fc+cfg.slope*(cfg.nSamples-1)/(2*cfg.fs);
fd=2*vr*fRef/cfg.c; delays=repelem((0:cfg.nTx-1)*cfg.slotTime,cfg.nRx);
positions=reshape((geom.rxOffset(:)+geom.txOffset),1,[]);
z=exp(-1i*2*pi*fRef/cfg.c*positions*sind(theta)).*exp(1i*2*pi*fd*delays);
zComp=z.*exp(-1i*2*pi*fd*delays);
angles=-60:0.05:60;
steer=exp(1i*2*pi*fRef/cfg.c*positions(:)*sind(angles));
raw=abs(z*steer); corrected=abs(zComp*steer);
[~,i]=max(raw); [~,j]=max(corrected);
demo=struct('angles',angles,'before',raw,'after',corrected,'trueAngle',theta, ...
    'estimatedBefore',angles(i),'estimatedAfter',angles(j),'fd',fd, ...
    'radialVelocity',vr,'tdmPhaseRad',2*pi*fd*cfg.slotTime);
end
