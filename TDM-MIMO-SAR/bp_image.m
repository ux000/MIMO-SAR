function [image,stat] = bp_image(rc,rangeAxis,geom,mimo,cfg,mode)
%BP_IMAGE ROI 内逐 chirp、全通道双站反投影。
% mode='correct'/'no_phase'/'trajectory_error'；三个分支使用同一含噪数据。
% BP 用真实各 Tx 时间戳，不再施加 MIMO 测角用的 TDM 对齐，避免重复补偿。
if nargin<6, mode='correct'; end
idx=find(mimo.mask); x=mimo.X(idx); y=mimo.Y(idx);
acc=complex(zeros(size(x))); sumMagnitude=zeros(size(x));
fRef=cfg.fc+cfg.slope*(cfg.nSamples-1)/(2*cfg.fs);
tStart=min(geom.time(:)); duration=max(geom.time(:))-tStart;
tAll=geom.time(:); apertureWeight=0.5-0.5*cos(2*pi*(tAll-tStart)/max(duration,eps));
apertureWeight=reshape(apertureWeight,size(geom.time));
totalWeight=0; nMeasurements=0; timer=tic;
for f=1:cfg.nFrames
    for m=1:cfg.nTx
        for l=1:cfg.nLoops
            t=geom.time(l,m,f); px=cfg.velocity(1)*t; py=cfg.velocity(2)*t;
            if strcmp(mode,'trajectory_error')
                py=py+cfg.trajectoryErrorM*sin(2*pi*3*(t-tStart)/max(duration,eps));
            end
            wt=apertureWeight(l,m,f);
            for r=1:cfg.nRx
                rt=hypot(x-px-geom.txOffset(m),y-py);
                rr=hypot(x-px-geom.rxOffset(r),y-py);
                tau=(rt+rr)/cfg.c;
                sample=interp1(rangeAxis,rc(:,l,m,r,f),cfg.c*tau/2,'linear',0);
                if ~strcmp(mode,'no_phase')
                    % 载频/快时间中心相位 + RVP 的共轭匹配。
                    sample=sample.*exp(-1i*2*pi*fRef*tau+1i*pi*cfg.slope*tau.^2);
                end
                acc=acc+wt*sample;
                sumMagnitude=sumMagnitude+wt*abs(sample);
                totalWeight=totalWeight+wt; nMeasurements=nMeasurements+1;
            end
        end
    end
end
image=complex(zeros(size(mimo.X))); image(idx)=acc/max(totalWeight,eps);
coherence=zeros(size(mimo.X)); coherence(idx)=abs(acc)./max(sumMagnitude,eps);
stat=struct('seconds',toc(timer),'measurements',nMeasurements,'pixels',numel(idx), ...
    'roiFraction',nnz(mimo.mask)/numel(mimo.mask),'coherence',coherence);
end
