function cfg = default_config()
%DEFAULT_CONFIG 参数集中设置；SI 单位。只需要基础 MATLAB，无额外工具箱。
cfg.c = 299792458;
cfg.fc = 77e9;                  % FMCW 起始频率，不直接采样 77 GHz 射频
cfg.slope = 21e12;              % 调频斜率，Hz/s
cfg.fs = 4e6;
cfg.nSamples = 64;
cfg.slotTime = 45e-6;           % 单个 Tx 的时间槽（含空闲时间）
cfg.nTx = 2; cfg.nRx = 4;
cfg.nLoops = 255;               % 一帧中每个 Tx 发射 255 次
cfg.nFrames = 13;
cfg.frameTime = 33.3e-3;
cfg.velocity = [1, 0];         % 已知平台速度，[沿轨 x，侧视 y]，m/s
cfg.snrDb = 10;                % 全 ADC 数据平均信号功率/平均 AWGN 功率
cfg.seed = 20261116;
cfg.rangeAttenuation = true;   % 幅度正比于 1/(Rtx*Rrx)
cfg.referenceRange = 5;
% 前两个点对应论文的 5 m、正负 0.5 度近邻点。
cfg.targets = [-5*sind(0.5),5*cosd(0.5); ...
                5*sind(0.5),5*cosd(0.5); -0.75,7.0; 0.85,8.2];
cfg.reflectivity = [1;1;0.9;0.9]; % 固定复散射系数；不要逐 chirp 随机相位
% 固定通道相位误差，顺序：Tx1-Rx1..4，然后 Tx2-Rx1..4。
cfg.channelPhaseDeg = [0,18,-22,30,14,-18,25,-12];
cfg.calibrateChannels = true;   % 使用已知标定值；并非从目标真值估计
cfg.nRangeFFT = 512;            % 零填充改善插值，不改善物理距离分辨率
cfg.nDopplerFFT = 256;
cfg.nAngleFFT = 256;
cfg.rangeLimits = [2,11];
cfg.angleLimitDeg = 60;
cfg.pfa = 1e-4;                % 近似模型下的 CA-CFAR 名义虚警概率
cfg.cfarTraining = [8,8];       % 距离、多普勒两维的训练环厚度
cfg.cfarGuard = [4,2];
cfg.maxDetections = 16;
cfg.detectionRelativeFloorDb = -23;
cfg.roiAngleWidthDeg = 10;      % 比论文 5 度更宽，给噪声和测角误差留余量
cfg.roiRangeWidth = 1.2;        % ROI 的 y 方向全宽，m
cfg.xLimits = [-1.6,1.6]; cfg.yLimits = [4.2,8.9];
cfg.dx = 0.01; cfg.dy = 0.025;
cfg.trajectoryErrorM = 1.2e-3;  % 对照组：y 方向正弦轨迹误差的峰值
cfg.figureVisible = 'on';
cfg.saveRaw = true;
cfg.saveFigures = true;
cfg.outputDir = fullfile(fileparts(mfilename('fullpath')),'results');
end
