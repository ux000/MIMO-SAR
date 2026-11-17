%% =========================================================================
% 77 GHz FMCW 同向 MIMO-SAR 墙面反射 NLOS 成像仿真
%
% 功能：
%   1. 汽车沿 x 方向运动，MIMO Tx/Rx 线阵也沿 x 方向排列（同向孔径）
%   2. 用 y=0 的有限反射墙构造 Tx->墙->隐藏目标->墙->Rx 双反射路径
%   3. x=4 的不透明遮挡墙阻断直达路径，反射路径从其下方绕行
%   4. 生成去斜后的 FMCW 多通道回波
%   5. 分别输出直达物理基准、直达模型失配诊断和墙面模型 NLOS-BP
%
% 重要建模说明：
%   - 这是一个二维 x-y 平面、静止点目标、理想镜面反射的基准仿真。
%   - y=0 是有限长度反射墙；x=4 的遮挡墙是另一段不透明障碍物。
%     这种“侧向反射面 + 拐角遮挡”的几何，才能同时满足：
%       (a) 目标直达路径被遮挡；
%       (b) 经 y=0 反射的路径可绕开遮挡物。
%   - Tx/Rx 实际坐标用于回波生成和成像，不用虚拟相位中心替代。
%   - 假设各 Tx 已理想分离（等效正交发射）；TDM 时应按实际发射时刻
%     更新平台位置并做运动补偿。
%   - 回波生成与 BP 使用完全相同的时延、去斜和相位约定。
%
% MATLAB：建议 R2019b 或更新版本
%% =========================================================================
clear; clc; close all;
rng(7);

%% 1. 雷达与平台参数
c       = 3e8;
fc      = 77e9;
lambda  = c/fc;

B       = 1.0e9;          % FMCW 带宽
Tchirp  = 40e-6;          % Chirp 时长
Kchirp  = B/Tchirp;       % 调频斜率
fs      = 10e6;           % ADC 采样率
Ns      = round(Tchirp*fs);
tFast   = (0:Ns-1)/fs;

Nfft    = 512;            % 距离频域 FFT 点数
fAxis   = (-Nfft/2:Nfft/2-1)*(fs/Nfft);
df      = fs/Nfft;

NT      = 4;
NR      = 4;
dTx     = lambda/2;
dRx     = lambda/2;
txOffset = ((0:NT-1)-(NT-1)/2)*dTx;  % Tx 沿 x 方向
rxOffset = ((0:NR-1)-(NR-1)/2)*dRx;  % Rx 沿 x 方向

Na      = 101;            % SAR 位置数
vPlatform = 5;            % 平台速度，m/s
Tr      = 4e-3;           % 相邻慢时间采样间隔
xStart  = -1.0;
xPlatform = xStart + (0:Na-1)*vPlatform*Tr;
yRadar  = 2.0;

SNR_dB  = 30;             % 接收复基带数据 SNR
Gamma   = 0.85*exp(1j*0.15); % 理想化墙面复反射系数
sigmaTarget = 1.0;        % 点目标复散射系数

fprintf('lambda = %.4f mm\n',lambda*1e3);
fprintf('SAR aperture = %.3f m\n',xPlatform(end)-xPlatform(1));
fprintf('MIMO Tx/Rx = %d x %d, channels = %d\n',NT,NR,NT*NR);

%% 2. 场景几何
% 反射墙：y=0，有限范围 x=[2,7]
yWall = 0;
wallXMin = 2.0;
wallXMax = 7.0;

% 不透明遮挡墙：x=4，y=[1,3]，阻断雷达到目标的直达线段
xBlock = 4.0;
blockYMin = 1.0;
blockYMax = 3.0;

% 隐藏目标（在遮挡墙右侧）
targetXY = [8.0, 2.0];

% 雷达轨迹 y=2；目标也位于 y=2
% 直达线经过 x=4,y约2，被遮挡墙挡住。
% 反射路径经 y=0，反射点在 x约3.5~4.5，低于遮挡墙。

%% 3. 检查目标对应的反射几何是否有效
xTxTest = xPlatform(ceil(Na/2)) + txOffset;
xRxTest = xPlatform(ceil(Na/2)) + rxOffset;
for m = 1:NT
    qTx = (xTxTest(m)+targetXY(1))/2;
    assert(qTx >= wallXMin && qTx <= wallXMax, ...
        'Tx 反射点落在有限反射墙范围外，请调整墙长或场景。');
end
for k = 1:NR
    qRx = (xRxTest(k)+targetXY(1))/2;
    assert(qRx >= wallXMin && qRx <= wallXMax, ...
        'Rx 反射点落在有限反射墙范围外，请调整墙长或场景。');
end

% 目标直达路径是否被遮挡（以中间平台位置、中间 Tx/Rx 检查）
xMid = xPlatform(ceil(Na/2));
xTxMid = xMid + txOffset(ceil(NT/2));
xRxMid = xMid + rxOffset(ceil(NR/2));
yAtBlockTx = yRadar + (targetXY(2)-yRadar)*(xBlock-xTxMid)/(targetXY(1)-xTxMid);
yAtBlockRx = yRadar + (targetXY(2)-yRadar)*(xBlock-xRxMid)/(targetXY(1)-xRxMid);
assert(yAtBlockTx >= blockYMin && yAtBlockTx <= blockYMax && ...
       yAtBlockRx >= blockYMin && yAtBlockRx <= blockYMax, ...
       '当前直达路径没有被遮挡墙挡住，请检查几何参数。');

%% 4. 绘制场景几何
figure('Color','w','Name','NLOS 场景几何');
hold on; grid on; axis equal;
plot([wallXMin wallXMax],[yWall yWall],'b-','LineWidth',5);
plot([xBlock xBlock],[blockYMin blockYMax],'k-','LineWidth',7);
plot(xPlatform,yRadar*ones(size(xPlatform)),'g.','MarkerSize',8);
plot(targetXY(1),targetXY(2),'rp','MarkerFaceColor','r','MarkerSize',14);

% 绘制中间通道的反射路径
m0 = ceil(NT/2); k0 = ceil(NR/2); n0 = ceil(Na/2);
xT0 = xPlatform(n0)+txOffset(m0);
xR0 = xPlatform(n0)+rxOffset(k0);
qT0 = [(xT0+targetXY(1))/2, yWall];
qR0 = [(xR0+targetXY(1))/2, yWall];
plot([xT0 qT0(1) targetXY(1) qR0(1) xR0], ...
     [yRadar yWall targetXY(2) yWall yRadar], ...
     'm--','LineWidth',1.6);
plot(qT0(1),qT0(2),'mo','MarkerFaceColor','m');
plot(qR0(1),qR0(2),'mo','MarkerFaceColor','m');

xlabel('x（沿轨方向，m）');
ylabel('y（横向方向，m）');
title('同向 MIMO-SAR 墙面反射 NLOS 场景');
legend('有限反射墙 y=0','不透明遮挡墙 x=4','雷达运动轨迹', ...
       '隐藏目标','NLOS 双反射路径','反射点','Location','best');
xlim([-1.5 9]); ylim([-0.5 3.5]);

%% 5. 生成 NLOS FMCW 去斜回波
% 数据维度：NT x NR x Na x Ns
rxIF = complex(zeros(NT,NR,Na,Ns,'single'));
tFastRow = single(tFast);

fprintf('正在生成 NLOS FMCW 回波...\n');
for m = 1:NT
    for k = 1:NR
        for n = 1:Na
            xTx = xPlatform(n)+txOffset(m);
            xRx = xPlatform(n)+rxOffset(k);
            xP  = targetXY(1);
            yP  = targetXY(2);

            % y=0 镜面反射：镜像法给出的反射点 x 坐标
            hMirror = yP + yRadar - 2*yWall;
            qTxX = xTx + (yRadar-yWall)/hMirror*(xP-xTx);
            qRxX = xRx + (yRadar-yWall)/hMirror*(xP-xRx);

            % 有限反射墙有效性检查
            if qTxX < wallXMin || qTxX > wallXMax || ...
               qRxX < wallXMin || qRxX > wallXMax
                error('存在通道的镜面反射点超出墙面范围。');
            end

            % Tx->wall->target 与 target->wall->Rx 的实际折线路径长度
            % 每一侧的两段长度相等（本例雷达与目标到墙的高度均为2m）
            % 镜像 Tx/Rx 到目标的直线距离，等于对应两段折线路径长度
            RtxLeg = sqrt((xP-xTx)^2 + hMirror^2);
            RrxLeg = sqrt((xP-xRx)^2 + hMirror^2);
            Lpath = RtxLeg+RrxLeg;
            tau = Lpath/c;

            % 传播幅度：简化球面扩散模型 + 两次墙面反射 + 目标散射
            alpha = sigmaTarget*(Gamma^2)/(max(RtxLeg*RrxLeg,1e-6));

            % 去斜信号：s_rx(t)*conj(s_tx(t))
            % beat 频率 f_b=-K*tau
            phaseFast = exp(-1j*2*pi*Kchirp*tau*tFastRow ...
                            +1j*pi*Kchirp*tau^2);
            echo = alpha*exp(-1j*2*pi*fc*tau).*phaseFast;

            rxIF(m,k,n,:) = single(echo);
        end
    end
end

% 按整体信号功率加入复高斯白噪声
sigPower = mean(abs(rxIF(:)).^2);
noisePower = sigPower/10^(SNR_dB/10);
noise = sqrt(noisePower/2)*(randn(size(rxIF),'single') + ...
                            1j*randn(size(rxIF),'single'));
rxIF_noiseOnly = noise; % 当前仿真没有直达目标回波，保存纯噪声供直达基准使用
rxIF = rxIF + noise;

fprintf('回波生成完成。平均信号功率 %.3e，噪声功率 %.3e\n', ...
        sigPower,noisePower);

%% 6. 距离频域处理
% 对快时间加窗并做 FFT。归一化使峰值量级不随 Ns 改变。
win = single(0.5-0.5*cos(2*pi*(0:Ns-1)/Ns));
win4 = reshape(win,[1 1 1 Ns]);
rangeFFT = fftshift(fft(rxIF.*win4,Nfft,4),4)/sum(win);
rangeFFTNoise = fftshift(fft(rxIF_noiseOnly.*win4,Nfft,4),4)/sum(win);

% 选取一个代表通道显示拍频谱
specShow = squeeze(abs(rangeFFT(ceil(NT/2),ceil(NR/2),ceil(Na/2),:)));
figure('Color','w','Name','FMCW 距离频谱');
plot(fAxis/1e6,20*log10(specShow/max(specShow)+eps),'LineWidth',1.2);
grid on; xlabel('拍频（MHz）'); ylabel('归一化幅度（dB）');
title('单通道去斜信号距离频谱');

%% 7. 成像网格
% 网格间隔取 2.5 cm，避免原 5 cm 网格过粗导致点扩散函数显示成块状
xGrid = 6.5:0.025:9.5;
yGrid = 0.8:0.025:3.2;
[XX,YY] = meshgrid(xGrid,yGrid);
[nY,nX] = size(XX);

imgNLOS = complex(zeros(nY,nX));          % 正确 NLOS 路径成像
imgDirect = complex(zeros(nY,nX));         % 纯噪声直达路径物理基准
imgMismatch = complex(zeros(nY,nX));       % NLOS 回波按直达模型处理的失配诊断

% 预先整理各通道的距离 FFT 数据，维度 Na x Nfft
RD = cell(NT,NR);
RDNoise = cell(NT,NR);
for m = 1:NT
    for k = 1:NR
        RD{m,k} = squeeze(rangeFFT(m,k,:,:));
        RDNoise{m,k} = squeeze(rangeFFTNoise(m,k,:,:));
    end
end

% 频率轴从负到正，便于线性插值
fMin = fAxis(1);

% SAR 慢时间加 Hann 孔径加权，降低有限孔径旁瓣
apWeight = 0.5-0.5*cos(2*pi*(0:Na-1)/(Na-1));
apWeight = apWeight / max(apWeight);

fprintf('正在执行直达模型失配成像与 NLOS-BP，请稍候...\n');

%% 8. BP 成像
for iy = 1:nY
    for ix = 1:nX
        xP = XX(iy,ix);
        yP = YY(iy,ix);

        sumN = 0;
        sumD = 0;
        sumMismatch = 0;
        nValidN = 0;
        nValidD = 0;
        nValidMismatch = 0;

        for m = 1:NT
            for k = 1:NR
                % 该通道在所有 SAR 位置的实际 Tx/Rx 坐标
                xTx = xPlatform + txOffset(m);
                xRx = xPlatform + rxOffset(k);

                %% 8.1 NLOS 路径：Tx->wall->pixel->wall->Rx
                hMirror = yP + yRadar - 2*yWall;
                qTxX = xTx + (yRadar-yWall)./hMirror.*(xP-xTx);
                qRxX = xRx + (yRadar-yWall)./hMirror.*(xP-xRx);

                validWall = qTxX>=wallXMin & qTxX<=wallXMax & ...
                            qRxX>=wallXMin & qRxX<=wallXMax;

                % 检查反射路径是否避开 x=xBlock 的不透明遮挡墙。
                % 路径的每一段与 x=xBlock 相交时，交点 y 必须在遮挡墙之外。
                validOcc = true(1,Na);
                for nn = 1:Na
                    % Tx -> qTx
                    if (xBlock-min(xTx(nn),qTxX(nn)) >= 0) && ...
                       (xBlock-max(xTx(nn),qTxX(nn)) <= 0)
                        if abs(qTxX(nn)-xTx(nn)) > eps
                            mu = (xBlock-xTx(nn))/(qTxX(nn)-xTx(nn));
                            yCross = yRadar + mu*(yWall-yRadar);
                            if yCross>=blockYMin && yCross<=blockYMax
                                validOcc(nn)=false;
                            end
                        end
                    end
                    % qTx -> pixel
                    if (xBlock-min(qTxX(nn),xP) >= 0) && ...
                       (xBlock-max(qTxX(nn),xP) <= 0)
                        if abs(xP-qTxX(nn)) > eps
                            mu = (xBlock-qTxX(nn))/(xP-qTxX(nn));
                            yCross = yWall + mu*(yP-yWall);
                            if yCross>=blockYMin && yCross<=blockYMax
                                validOcc(nn)=false;
                            end
                        end
                    end
                    % pixel -> qRx
                    if (xBlock-min(xP,qRxX(nn)) >= 0) && ...
                       (xBlock-max(xP,qRxX(nn)) <= 0)
                        if abs(qRxX(nn)-xP) > eps
                            mu = (xBlock-xP)/(qRxX(nn)-xP);
                            yCross = yP + mu*(yWall-yP);
                            if yCross>=blockYMin && yCross<=blockYMax
                                validOcc(nn)=false;
                            end
                        end
                    end
                    % qRx -> Rx
                    if (xBlock-min(qRxX(nn),xRx(nn)) >= 0) && ...
                       (xBlock-max(qRxX(nn),xRx(nn)) <= 0)
                        if abs(xRx(nn)-qRxX(nn)) > eps
                            mu = (xBlock-qRxX(nn))/(xRx(nn)-qRxX(nn));
                            yCross = yWall + mu*(yRadar-yWall);
                            if yCross>=blockYMin && yCross<=blockYMax
                                validOcc(nn)=false;
                            end
                        end
                    end
                end

                % 镜像源到像素的距离，等于 Tx->墙->像素 / 像素->墙->Rx
                Rtx = sqrt((xP-xTx).^2+hMirror.^2);
                Rrx = sqrt((xP-xRx).^2+hMirror.^2);
                tauN = (Rtx+Rrx)/c;

                % NLOS 路径的拍频查询位置
                fQueryN = -Kchirp*tauN;
                valN = interp_range_fft(RD{m,k},fQueryN,fMin,df,Nfft);

                % 与回波生成完全一致的去斜相位补偿
                phaseCompN = exp(1j*2*pi*fc*tauN - 1j*pi*Kchirp*tauN.^2);
                useN = validWall & validOcc;
                sumN = sumN + sum(valN(useN).*phaseCompN(useN).*apWeight(useN));
                nValidN = nValidN + sum(apWeight(useN));

                %% 8.2 直达路径物理基准 + NLOS失配诊断（严格分开）
                RtxD = sqrt((xP-xTx).^2+(yP-yRadar)^2);
                RrxD = sqrt((xP-xRx).^2+(yP-yRadar)^2);
                tauD = (RtxD+RrxD)/c;
                fQueryD = -Kchirp*tauD;
                phaseCompD = exp(1j*2*pi*fc*tauD ...
                                 -1j*pi*Kchirp*tauD.^2);

                % (A) 直达路径可见性：Tx->像素和像素->Rx均不能穿过遮挡墙
                visibleTx = segment_clear_vertical_blocker( ...
                    xTx,yRadar,xP,yP,xBlock,blockYMin,blockYMax);
                visibleRx = segment_clear_vertical_blocker( ...
                    xP,yP,xRx,yRadar,xBlock,blockYMin,blockYMax);
                useD = visibleTx & visibleRx;

                % 当前场景只生成 NLOS 目标回波，没有直达目标回波。
                % 因此物理直达基准必须使用纯噪声数据。
                valD_noise = interp_range_fft( ...
                    RDNoise{m,k},fQueryD,fMin,df,Nfft);
                sumD = sumD + sum(valD_noise(useD).* ...
                                  phaseCompD(useD).*apWeight(useD));
                nValidD = nValidD + sum(apWeight(useD));

                % (B) 独立的失配诊断：将含 NLOS 回波的数据按直达模型处理，
                % 不加可见性掩膜，仅用于显示错误传播模型造成的伪影。
                valMismatch = interp_range_fft( ...
                    RD{m,k},fQueryD,fMin,df,Nfft);
                sumMismatch = sumMismatch + sum( ...
                    valMismatch.*phaseCompD.*apWeight);
                nValidMismatch = nValidMismatch + sum(apWeight);
            end
        end

        if nValidN>0
            imgNLOS(iy,ix)=sumN/nValidN;
        end
        if nValidD>0
            imgDirect(iy,ix)=sumD/nValidD;
        else
            imgDirect(iy,ix)=0; % 所有直达观测均被遮挡
        end
        if nValidMismatch>0
            imgMismatch(iy,ix)=sumMismatch/nValidMismatch;
        end
    end
    if mod(iy,5)==0
        fprintf('成像进度：%d / %d 行\n',iy,nY);
    end
end

%% 9. 归一化与显示
% 统一以 NLOS-BP 峰值作为参考，避免每幅图单独归一化造成误导。
commonScale = max(abs(imgNLOS(:))) + eps;
imgN_dB = 20*log10(abs(imgNLOS)/commonScale+1e-6);
imgD_dB = 20*log10(abs(imgDirect)/commonScale+1e-6);
imgMismatch_dB = 20*log10(abs(imgMismatch)/commonScale+1e-6);

% 图A：物理直达路径基准。当前仿真无直达目标回波，因此使用纯噪声数据。
figure('Color','w','Name','直达路径物理基准（纯噪声）');
imagesc(xGrid,yGrid,imgD_dB);
axis xy equal tight; colorbar; caxis([-60 0]);
hold on; plot(targetXY(1),targetXY(2),'rp', ...
    'MarkerFaceColor','r','MarkerSize',12);
xlabel('x（m）'); ylabel('y（m）');
title('直达路径基准：纯噪声数据 + 直达可见性掩膜');

% 图B：失配诊断。NLOS回波按直达路径延迟/相位处理，不施加可见性掩膜。
figure('Color','w','Name','直达模型失配诊断');
imagesc(xGrid,yGrid,imgMismatch_dB);
axis xy equal tight; colorbar; caxis([-40 0]);
hold on; plot(targetXY(1),targetXY(2),'rp', ...
    'MarkerFaceColor','r','MarkerSize',12);
xlabel('x（m）'); ylabel('y（m）');
title('NLOS回波按直达模型处理（失配伪影，非真实直达成像）');

% 图C：正确的墙面反射 NLOS-BP
figure('Color','w','Name','墙面反射 NLOS-BP');
imagesc(xGrid,yGrid,imgN_dB);
axis xy equal tight; colorbar; caxis([-40 0]);
hold on; plot(targetXY(1),targetXY(2),'rp', ...
    'MarkerFaceColor','r','MarkerSize',12);
xlabel('x（m）'); ylabel('y（m）');
title('墙面反射 NLOS-BP 成像结果');

% 图D：三种结果并列，统一幅度参考
figure('Color','w','Name','三种处理结果对比');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile;
imagesc(xGrid,yGrid,imgD_dB); axis xy equal tight; colorbar; caxis([-60 0]);
hold on; plot(targetXY(1),targetXY(2),'rp','MarkerFaceColor','r');
xlabel('x（m）'); ylabel('y（m）');
title('直达基准（纯噪声）');

nexttile;
imagesc(xGrid,yGrid,imgMismatch_dB); axis xy equal tight; colorbar; caxis([-40 0]);
hold on; plot(targetXY(1),targetXY(2),'rp','MarkerFaceColor','r');
xlabel('x（m）'); ylabel('y（m）');
title('直达模型失配诊断');

nexttile;
imagesc(xGrid,yGrid,imgN_dB); axis xy equal tight; colorbar; caxis([-40 0]);
hold on; plot(targetXY(1),targetXY(2),'rp','MarkerFaceColor','r');
xlabel('x（m）'); ylabel('y（m）');
title('NLOS-BP');

sgtitle('同向 MIMO-SAR 非视域成像：直达基准、失配诊断与 NLOS-BP');

%% 10. 输出峰值定位误差
[~,idxMax] = max(abs(imgNLOS(:)));
[iyMax,ixMax] = ind2sub(size(imgNLOS),idxMax);
peakXY = [xGrid(ixMax),yGrid(iyMax)];
locError = norm(peakXY-targetXY);

fprintf('\n========== 仿真结果 ==========\n');
fprintf('真实目标位置：(%.3f, %.3f) m\n',targetXY(1),targetXY(2));
fprintf('NLOS-BP 峰值位置：(%.3f, %.3f) m\n',peakXY(1),peakXY(2));
fprintf('定位误差：%.4f m\n',locError);
fprintf('NLOS-BP 峰值幅度（归一化前）：%.4e\n',abs(imgNLOS(iyMax,ixMax)));
fprintf('直达物理基准峰值 / NLOS-BP 峰值：%.2f dB\n', ...
    20*log10(max(abs(imgDirect(:))) / (max(abs(imgNLOS(:)))+eps) + eps));
fprintf('直达模型失配峰值 / NLOS-BP 峰值：%.2f dB\n', ...
    20*log10(max(abs(imgMismatch(:))) / (max(abs(imgNLOS(:)))+eps) + eps));
fprintf('==============================\n');

%% 11. 说明
disp('仿真完成：直达物理基准使用纯噪声数据；失配诊断单独使用 NLOS 回波。');
disp('注意：本程序是理想二维镜面反射基准模型，不包含墙面粗糙散射、');
disp('多径杂波、目标运动、天线方向图、互耦和TDM时序误差。');

%% ======================== 局部函数 =========================

function isClear = segment_clear_vertical_blocker(x1,y1,x2,y2,xWall,yMin,yMax)
% 判断二维线段 (x1,y1)->(x2,y2) 是否穿过有限竖直遮挡墙：
%   x=xWall, y in [yMin,yMax]。
% 输入坐标可以是标量或等长向量，输出为对应逻辑行向量。
    x1 = x1(:).'; y1 = y1(:).';
    x2 = x2(:).'; y2 = y2(:).';

    n = max([numel(x1),numel(y1),numel(x2),numel(y2)]);
    if ~(numel(x1)==1 || numel(x1)==n) || ...
       ~(numel(y1)==1 || numel(y1)==n) || ...
       ~(numel(x2)==1 || numel(x2)==n) || ...
       ~(numel(y2)==1 || numel(y2)==n)
        error('遮挡判断函数的输入必须为标量或相同长度的向量。');
    end

    if isscalar(x1), x1 = repmat(x1,1,n); end
    if isscalar(y1), y1 = repmat(y1,1,n); end
    if isscalar(x2), x2 = repmat(x2,1,n); end
    if isscalar(y2), y2 = repmat(y2,1,n); end

    dx = x2-x1;
    crossesX = (xWall >= min(x1,x2)) & (xWall <= max(x1,x2));
    yCross = nan(1,n);
    nz = abs(dx)>eps;
    yCross(nz) = y1(nz)+(y2(nz)-y1(nz)).*(xWall-x1(nz))./dx(nz);

    % 若线段与遮挡墙所在 x 坐标共线，且 y 范围有重叠，也视为遮挡
    collinear = ~nz & abs(x1-xWall)<=eps & ...
        max(min(y1,y2),yMin)<=min(max(y1,y2),yMax);

    blocked = (crossesX & nz & yCross>=yMin & yCross<=yMax) | collinear;
    isClear = ~blocked;
end

function yq = interp_range_fft(RDmat,fQuery,fMin,df,Nfft)
% 对每个慢时间位置，在距离 FFT 频率轴上进行复数线性插值。
% RDmat: Na x Nfft；fQuery: 1 x Na 或 Na x 1
    fQuery = fQuery(:);
    NaLocal = size(RDmat,1);
    if numel(fQuery) ~= NaLocal
        error('频率查询向量长度与慢时间采样数不一致。');
    end

    u = (fQuery-fMin)/df + 1;
    i0 = floor(u);
    frac = u-i0;

    valid = i0>=1 & i0<Nfft & isfinite(u);
    yq = complex(zeros(NaLocal,1));

    rows = find(valid);
    if ~isempty(rows)
        cols0 = i0(rows);
        cols1 = cols0+1;
        ind0 = sub2ind(size(RDmat),rows,cols0);
        ind1 = sub2ind(size(RDmat),rows,cols1);
        a = RDmat(ind0);
        b = RDmat(ind1);
        yq(rows) = a.*(1-frac(rows))+b.*frac(rows);
    end
    yq = yq.';
end
