function mimo = mimo_roi(rc,rangeAxis,geom,cfg)
%MIMO_ROI 第一帧：距离/多普勒FFT -> CA-CFAR -> 峰聚类 -> TDM补偿 -> 角度FFT。
% 本函数不读取 cfg.targets。ROI 完全由有噪声的观测生成。
assert(cfg.nDopplerFFT>=cfg.nLoops,'多普勒FFT点数不足。');
n = (0:cfg.nLoops-1)';
wd = 0.5-0.5*cos(2*pi*n/(cfg.nLoops-1));
v = fftshift(fft(rc(:,:,:,:,1).*reshape(wd,1,[]),cfg.nDopplerFFT,2),2)/sum(wd);
powerRV = sum(sum(abs(v).^2,4),3);
fd = (-floor(cfg.nDopplerFFT/2):ceil(cfg.nDopplerFFT/2)-1)/(cfg.nDopplerFFT*geom.loopTime);
% 功率非相干叠加；独立Gamma噪声下的F比例阈值。窗/零填充会引入相关，
% 因而 pfa 是名义值，不能将它当作实际测得的虚警概率。
tr=cfg.cfarTraining; gd=cfg.cfarGuard; outer=tr+gd;
kernel=ones(2*outer(1)+1,2*outer(2)+1);
kernel(outer(1)-gd(1)+1:outer(1)+gd(1)+1,outer(2)-gd(2)+1:outer(2)+gd(2)+1)=0;
nTrain=sum(kernel(:)); channels=cfg.nTx*cfg.nRx;
beta=betaincinv(1-cfg.pfa,channels,channels*nTrain);
alpha=nTrain*beta/(1-beta);
noiseEstimate=conv2(powerRV,kernel,'same')/nTrain;
hit=powerRV>alpha*noiseEstimate;
hit(1:outer(1),:)=false; hit(end-outer(1)+1:end,:)=false;
hit(:,1:outer(2))=false; hit(:,end-outer(2)+1:end)=false;
hit=hit & powerRV>max(powerRV(:))*10^(cfg.detectionRelativeFloorDb/10);
% 3x3 局部峰分组；无 Image Processing Toolbox 依赖。
for a=-1:1
    for b=-1:1
        if a~=0 || b~=0, hit=hit & powerRV>=circshift(powerRV,[a,b]); end
    end
end
[ir,id]=find(hit); [~,order]=sort(powerRV(hit),'descend');
order=order(1:min(cfg.maxDetections,numel(order))); ir=ir(order); id=id(order);
fRef=cfg.fc+cfg.slope*(cfg.nSamples-1)/(2*cfg.fs);
spatialFreq=(-floor(cfg.nAngleFFT/2):ceil(cfg.nAngleFFT/2)-1)/cfg.nAngleFFT;
sinAngle=-spatialFreq*cfg.c/(fRef*(cfg.c/cfg.fc/2));
valid=abs(sinAngle)<=sind(cfg.angleLimitDeg);
angles=asind(sinAngle(valid)); [angles,sortAngle]=sort(angles);
% 每个 Tx 块先补偿同一 Doppler bin 的异步发射相位，再构造虚拟阵列。
anglePower=zeros(numel(rangeAxis),numel(angles));
for d=1:numel(fd)
    vv=reshape(permute(v(:,d,:,:),[1,4,3,2]),numel(rangeAxis),[]);
    delays=repelem((0:cfg.nTx-1)*cfg.slotTime,cfg.nRx);
    vv=vv.*exp(-1i*2*pi*fd(d)*delays);
    aa=fftshift(fft(vv,cfg.nAngleFFT,2),2);
    aa=aa(:,valid); aa=aa(:,sortAngle);
    anglePower=anglePower+abs(aa).^2;
end
tRef=mean(geom.time(:,1,1)); radarRef=cfg.velocity*tRef;
det=zeros(numel(ir),7);
for k=1:numel(ir)
    vv=reshape(permute(v(ir(k),id(k),:,:),[4,3,1,2]),1,[]);
    vv=vv.*exp(-1i*2*pi*fd(id(k))*delays);
    aa=fftshift(fft(vv,cfg.nAngleFFT)); aa=aa(valid); aa=aa(sortAngle);
    [~,ia]=max(abs(aa)); theta=angles(ia); range=rangeAxis(ir(k));
    det(k,:)=[range,fd(id(k))*cfg.c/(2*fRef),theta, ...
        radarRef(1)+range*sind(theta),radarRef(2)+range*cosd(theta), ...
        powerRV(ir(k),id(k)),fd(id(k))];
end
x=cfg.xLimits(1):cfg.dx:cfg.xLimits(2); y=cfg.yLimits(1):cfg.dy:cfg.yLimits(2);
[X,Y]=meshgrid(x,y); mask=false(size(X)); boxes=zeros(size(det,1),4);
for k=1:size(det,1)
    halfX=max(0.15,det(k,1)*deg2rad(cfg.roiAngleWidthDeg)/2);
    halfY=cfg.roiRangeWidth/2;
    boxes(k,:)=[det(k,4)-halfX,det(k,5)-halfY,2*halfX,2*halfY];
    mask=mask | (abs(X-det(k,4))<=halfX & abs(Y-det(k,5))<=halfY);
end
if ~any(mask(:))
    error('MIMO未检出可用ROI：检查信噪比、CFAR参数或场景范围；程序不会用目标真值补ROI。');
end
mimo=struct('rangeAxis',rangeAxis,'dopplerHz',fd,'powerRV',powerRV, ...
    'anglesDeg',angles,'anglePower',anglePower,'detections',det, ...
    'x',x,'y',y,'X',X,'Y',Y,'mask',mask,'boxes',boxes, ...
    'cfarAlpha',alpha,'referenceTime',tRef);
% detections各列：[距离m，远离为正径向速度m/s，角度deg，全局x，全局y，功率，fD Hz]
end
