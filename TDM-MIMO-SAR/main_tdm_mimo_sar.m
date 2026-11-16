function out = main_tdm_mimo_sar(overrides)
%MAIN_TDM_MIMO_SAR 带噪声、中文图的 TDM-MIMO-SAR 完整教学仿真。
% 直接运行：out = main_tdm_mimo_sar;
% 修改参数：out = main_tdm_mimo_sar(struct('snrDb',0));
% 采用论文的 MIMO->ROI->BP 空间层级思路。此版本使用逐chirp精确双站BP，
% 不宣称复现论文 K-chirp RVA 压缩加速或雷达里程计/RANSAC。
cfg=default_config();
if nargin>0
    names=fieldnames(overrides);
    for k=1:numel(names)
        assert(isfield(cfg,names{k}),'未知参数：%s',names{k});
        cfg.(names{k})=overrides.(names{k});
    end
end
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
fprintf('1/5 生成真实时序的 TDM 双站回波，加入 AWGN ...\n');
[adc,geom,echoInfo]=simulate_echo(cfg);
fprintf('    实测 SNR %.3f dB，距离分辨率 %.3f m，合成孔径 %.3f m\n', ...
    echoInfo.measuredSnrDb,echoInfo.rangeResolution,echoInfo.apertureLength);
[rc,rangeAxis,rangeMeta]=range_compress(adc,cfg);
fprintf('2/5 MIMO 检测、TDM相位对齐、角度FFT及ROI ...\n');
mimo=mimo_roi(rc,rangeAxis,geom,cfg);
fprintf('    检出 %d 个峰，ROI占全图 %.1f%%\n',size(mimo.detections,1),100*nnz(mimo.mask)/numel(mimo.mask));
fprintf('3/5 ROI 内做逐chirp精确BP，及两个相位对照组 ...\n');
[imCorrect,bpStat]=bp_image(rc,rangeAxis,geom,mimo,cfg,'correct');
[imNoPhase,noPhaseStat]=bp_image(rc,rangeAxis,geom,mimo,cfg,'no_phase');
[imTrajectoryError,errorStat]=bp_image(rc,rangeAxis,geom,mimo,cfg,'trajectory_error');
demo=tdm_phase_demo(cfg,geom);
out=struct('cfg',cfg,'geom',geom,'echoInfo',echoInfo,'rangeMeta',rangeMeta, ...
    'rangeAxis',rangeAxis,'mimo',mimo,'image',imCorrect,'imageNoPhase',imNoPhase, ...
    'imageTrajectoryError',imTrajectoryError,'bpStat',bpStat, ...
    'noPhaseStat',noPhaseStat,'errorStat',errorStat,'tdmDemo',demo);
fprintf('4/5 生成中文图及数值报告 ...\n');
make_plots(adc,rc,out);
write_report(out);
fprintf('5/5 保存结果 ...\n');
save(fullfile(cfg.outputDir,'simulation_results.mat'),'out','-v7.3');
if cfg.saveRaw, save(fullfile(cfg.outputDir,'echo_data.mat'),'adc','geom','cfg','-v7.3'); end
fprintf('完成。结果路径：%s\n',cfg.outputDir);
end

function write_report(out)
cfg=out.cfg; f=fopen(fullfile(cfg.outputDir,'数值结果.txt'),'w','n','UTF-8');
cleanup=onCleanup(@()fclose(f));
fprintf(f,'TDM-MIMO-SAR 数值结果\n运行环境：MATLAB %s\n',version);
fprintf(f,'设定/实测SNR：%.3f / %.3f dB\n',cfg.snrDb,out.echoInfo.measuredSnrDb);
fprintf(f,'带宽 %.3f MHz，物理距离分辨率 %.4f m\n',out.echoInfo.bandwidth/1e6,out.echoInfo.rangeResolution);
fprintf(f,'合成孔径 %.4f m；单槽平台位移 %.6f m\n',out.echoInfo.apertureLength,out.echoInfo.maxSingleSlotDisplacement);
fprintf(f,'ROI像素 %d / %d (%.2f%%)；仅像素数量减少 %.2f 倍\n',nnz(out.mimo.mask),numel(out.mimo.mask),100*out.bpStat.roiFraction,1/out.bpStat.roiFraction);
fprintf(f,'该比值不是总运行时间加速比；本程序没有执行论文K-chirp压缩。\n');
fprintf(f,'三次BP用时（正确/无补偿/轨迹误差）：%.2f / %.2f / %.2f 秒\n',out.bpStat.seconds,out.noPhaseStat.seconds,out.errorStat.seconds);
fprintf(f,'TDM受控单目标角：真值 %.2f，补偿前 %.2f，补偿后 %.2f 度；跨Tx运动相位 %.3f rad\n',out.tdmDemo.trueAngle,out.tdmDemo.estimatedBefore,out.tdmDemo.estimatedAfter,out.tdmDemo.tdmPhaseRad);
fprintf(f,'目标真值仅用于以下验收与叠图，不用于ROI生成：\n');
for k=1:size(cfg.targets,1)
    [~,ix]=min(abs(out.mimo.x-cfg.targets(k,1))); [~,iy]=min(abs(out.mimo.y-cfg.targets(k,2)));
    fprintf(f,'目标%d (%.4f,%.4f)m，真值栅格是否在ROI：%d，相干系数 %.4f\n',k,cfg.targets(k,1),cfg.targets(k,2),out.mimo.mask(iy,ix),out.bpStat.coherence(iy,ix));
end
fprintf(f,'检测列：距离m / 远离为正速度m/s / 角度deg / 全局xm / 全局ym\n');
fprintf(f,'%10.4f %10.4f %10.4f %10.4f %10.4f\n',out.mimo.detections(:,1:5).');
end
