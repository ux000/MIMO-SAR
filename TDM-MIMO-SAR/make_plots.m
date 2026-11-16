function make_plots(adc,rc,out)
%MAKE_PLOTS 所有坐标、图题和色条采用中文；图中真值只用于显示与验收。
cfg=out.cfg; mi=out.mimo; fdir=cfg.outputDir;
font='Microsoft YaHei';
available=listfonts;
if ~any(strcmpi(available,font))
    if any(strcmpi(available,'SimHei')), font='SimHei'; else, font=get(groot,'defaultAxesFontName'); end
end
f=newfigure('回波与距离压缩',cfg);
tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
nexttile;
slots=0:6; stairs(slots*cfg.slotTime*1e6,mod(slots,cfg.nTx)+1,'LineWidth',1.6);
ylim([0.7,cfg.nTx+0.3]); yticks(1:cfg.nTx); grid on;
xlabel('慢时间（微秒）'); ylabel('发射天线编号'); title('TDM 时序：每个时间槽仅一个 Tx 发射');
nexttile;
tf=(0:cfg.nSamples-1)'/cfg.fs;
plot(tf*1e6,real(adc(:,1,1,1,1)),'LineWidth',1); hold on;
plot(tf*1e6,imag(adc(:,1,1,1,1)),'LineWidth',1); grid on;
xlabel('快时间（微秒）'); ylabel('去斜回波幅度'); title(sprintf('首个 chirp 的含噪复回波，SNR = %.0f dB',cfg.snrDb)); legend('同相分量','正交分量');
nexttile;
plot(out.rangeAxis,20*log10(max(abs(rc(:,1,1,1,1)),1e-9)),'LineWidth',1.4); grid on;
xlabel('距离（米）'); ylabel('距离压缩幅度（分贝）'); title(sprintf('加窗距离 FFT，物理分辨率 %.3f 米',out.echoInfo.rangeResolution));
nexttile;
imagesc((0:cfg.nLoops-1)*out.geom.loopTime*1e3,out.rangeAxis,relative_db(abs(rc(:,:,1,1,1)))); axis xy;
xlabel('帧内慢时间（毫秒）'); ylabel('距离（米）'); title('首帧单通道的距离—慢时间回波'); caxis([-40,0]); chinese_colorbar('相对幅度（分贝）');
finish(f,font,fullfile(fdir,'01_回波与距离压缩.png'),cfg);

f=newfigure('MIMO检测与ROI',cfg); tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
nexttile;
fRef=out.rangeMeta.fRef;
imagesc(mi.dopplerHz*cfg.c/(2*fRef),out.rangeAxis,relative_db(sqrt(mi.powerRV))); axis xy;
xlabel('径向速度（米/秒，远离为正）'); ylabel('距离（米）'); title('多通道非相干距离—多普勒图'); caxis([-40,0]); chinese_colorbar('相对功率（分贝）');
hold on; plot(mi.detections(:,2),mi.detections(:,1),'wo','MarkerSize',7);
nexttile;
imagesc(mi.anglesDeg,out.rangeAxis,relative_db(sqrt(mi.anglePower))); axis xy;
xlabel('方位角（度）'); ylabel('距离（米）'); title('TDM 补偿后的 MIMO 距离—角度图'); caxis([-35,0]); chinese_colorbar('相对功率（分贝）');
hold on; plot(mi.detections(:,3),mi.detections(:,1),'wo','MarkerSize',6);
nexttile;
imagesc(mi.x,mi.y,double(mi.mask)); axis xy equal tight; colormap(gca,[0.95,0.96,0.98;0.72,0.85,0.98]); hold on;
for k=1:size(mi.boxes,1), rectangle('Position',mi.boxes(k,:),'EdgeColor',[0.15,0.4,0.75],'LineWidth',0.9); end
plot(mi.detections(:,4),mi.detections(:,5),'b+','LineWidth',1.2);
plot(cfg.targets(:,1),cfg.targets(:,2),'rx','LineWidth',1.5,'MarkerSize',9);
xlabel('沿轨位置（米）'); ylabel('侧视距离（米）'); title(sprintf('MIMO 自动 ROI，占全图 %.1f%%',100*out.bpStat.roiFraction));
finish(f,font,fullfile(fdir,'02_MIMO检测与ROI.png'),cfg);

f=newfigure('BP相位补偿对比',cfg); tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
reference=max(abs(out.image(:)));
ims={out.imageNoPhase,out.imageTrajectoryError,out.image};
names={'未做 SAR 相位补偿',sprintf('轨迹存在 %.1f 毫米正弦误差',cfg.trajectoryErrorM*1e3),'正确轨迹、载频及 RVP 补偿'};
for k=1:3
    nexttile;
    db=20*log10(max(abs(ims{k})/max(reference,eps),1e-8)); db(~mi.mask)=NaN;
    h=imagesc(mi.x,mi.y,db); set(h,'AlphaData',mi.mask); set(gca,'Color',[0.93,0.94,0.95]); axis xy equal tight;
    hold on; plot(cfg.targets(:,1),cfg.targets(:,2),'wx','LineWidth',1.1);
    xlabel('沿轨位置（米）'); ylabel('侧视距离（米）'); title(names{k}); caxis([-40,0]);
    chinese_colorbar('相对幅度（分贝，共同参考）');
end
finish(f,font,fullfile(fdir,'03_BP相位补偿对比.png'),cfg);

f=newfigure('近邻目标剖面',cfg); tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
[~,iy]=min(abs(mi.y-cfg.targets(1,2)));
sel=abs(mi.x)<0.25;
nexttile;
plot(mi.x(sel),relative_db(abs(out.image(iy,sel))),'LineWidth',1.8); hold on;
[~,ir]=min(abs(out.rangeAxis-hypot(cfg.targets(1,1),cfg.targets(1,2))));
xAngle=cfg.velocity(1)*mi.referenceTime+out.rangeAxis(ir)*sind(mi.anglesDeg);
ra=sqrt(mi.anglePower(ir,:)); ra=ra/max(ra);
plot(xAngle,20*log10(max(ra,1e-8)),'--','LineWidth',1.5);
for k=1:min(2,size(cfg.targets,1)), xline(cfg.targets(k,1),':','Color',[0.4,0.4,0.4]); end
xlim([-0.25,0.25]); ylim([-40,2]); grid on;
xlabel('沿轨位置（米）'); ylabel('各自归一化幅度（分贝）'); title('近邻目标的横向剖面'); legend('ROI 内相干 BP','单帧 MIMO','Location','southwest');
nexttile;
db=relative_db(abs(out.image)); db(~mi.mask)=NaN;
imagesc(mi.x,mi.y,db); axis xy; xlim([-0.25,0.25]); ylim([4.5,5.5]); caxis([-35,0]);
hold on; plot(cfg.targets(1:min(2,end),1),cfg.targets(1:min(2,end),2),'wx','LineWidth',1.4);
xlabel('沿轨位置（米）'); ylabel('侧视距离（米）'); title('约 8.73 厘米间隔的两个静止点'); chinese_colorbar('相对幅度（分贝）');
finish(f,font,fullfile(fdir,'04_近邻目标剖面.png'),cfg);

f=newfigure('TDM相位补偿',cfg); tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact'); demo=out.tdmDemo;
nexttile;
plot(demo.angles,20*log10(max(demo.before/(cfg.nTx*cfg.nRx),1e-8)),'--','LineWidth',1.4); hold on;
plot(demo.angles,20*log10(max(demo.after/(cfg.nTx*cfg.nRx),1e-8)),'LineWidth',1.8);
xline(demo.trueAngle,':','真实角度'); xlim([-5,45]); ylim([-35,1]); grid on;
xlabel('方位角（度）'); ylabel('波束幅度（分贝）'); title('独立单目标受控测角实验'); legend('未补偿 TDM','已补偿 TDM','Location','southwest');
nexttile;
ch=1:cfg.nTx*cfg.nRx; residual=2*pi*demo.fd*repelem((0:cfg.nTx-1)*cfg.slotTime,cfg.nRx);
stem(ch,residual,'filled','LineWidth',1.4); hold on; plot(ch,zeros(size(ch)),'o-','LineWidth',1.2); grid on;
xlabel('虚拟通道编号'); ylabel('运动引入的相位（弧度）'); title(sprintf('径向速度 %.1f 米/秒，跨 Tx 相位 %.2f 弧度',demo.radialVelocity,demo.tdmPhaseRad));
legend('补偿前','补偿后');
finish(f,font,fullfile(fdir,'05_TDM相位补偿.png'),cfg);

wave=fmcw_waveform_demo(cfg,out.geom);
f=newfigure('发射接收与混频',cfg); tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
select=wave.time<=0.4e-6; t=wave.time(select); tau=wave.tau;
nexttile;
plot(t*1e6,cfg.slope*t/1e6,'LineWidth',1.5); hold on;
valid=t>=tau; plot(t(valid)*1e6,cfg.slope*(t(valid)-tau)/1e6,'--','LineWidth',1.5); grid on;
xlabel('快时间（微秒）'); ylabel('相对载频的瞬时频率（兆赫）'); title(sprintf('首个目标的双程延迟 %.2f 纳秒',tau*1e9)); legend('发射 chirp','延迟接收 chirp','Location','northwest');
nexttile;
plot(t*1e6,real(wave.txBaseband(select)),'LineWidth',1.2); hold on;
plot(t*1e6,real(wave.rxBaseband(select)),'--','LineWidth',1.2); grid on;
xlabel('快时间（微秒）'); ylabel('复包络实部（归一化）'); title('高采样率显式构造发射与接收波形'); legend('发射复包络','接收复包络','Location','southwest');
nexttile;
select=wave.time>=tau & wave.time<=8e-6;
plot(wave.time(select)*1e6,real(wave.ifSignal(select)),'LineWidth',1.2); grid on;
xlabel('快时间（微秒）'); ylabel('去斜信号实部（归一化）'); title(sprintf('发射乘接收共轭，拍频 %.3f 兆赫',wave.beatHz/1e6));
finish(f,font,fullfile(fdir,'06_发射接收与混频.png'),cfg);
plot_scene_geometry(cfg,out.geom);
end

function f=newfigure(name,cfg)
f=figure('Name',name,'Color','w','Position',[60,80,1500,620],'Visible',cfg.figureVisible);
end
function db=relative_db(a)
db=20*log10(max(a/max(max(a(:)),eps),1e-8));
end
function chinese_colorbar(label)
cb=colorbar; cb.Label.String=label;
end
function finish(f,font,path,cfg)
set(findall(f,'-property','FontName'),'FontName',font);
set(findall(f,'-property','FontSize'),'FontSize',11);
if cfg.saveFigures, exportgraphics(f,path,'Resolution',170); end
if strcmpi(cfg.figureVisible,'off'), close(f); end
end
