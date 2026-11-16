function [rc, rangeAxis, meta] = range_compress(adc,cfg)
%RANGE_COMPRESS 加窗、距离FFT、窗中心相位校正和固定通道标定。
% exp(+j*2*pi*f*tCenter) 使复数插值不会跨越未校正的窗中心相位斜坡。
% 经此中心化后，有效载频为 fc+slope*tCenter，BP 必须用同一约定。
assert(cfg.nRangeFFT>=cfg.nSamples,'距离FFT点数必须大于等于采样点数。');
n = (0:cfg.nSamples-1)';
w = 0.5-0.5*cos(2*pi*n/(cfg.nSamples-1));
freq = (0:cfg.nRangeFFT-1)'*cfg.fs/cfg.nRangeFFT;
meta.tCenter = (cfg.nSamples-1)/(2*cfg.fs);
meta.fRef = cfg.fc+cfg.slope*meta.tCenter;
meta.windowGain = sum(w);
allRange = freq*cfg.c/(2*cfg.slope);
keep = allRange>=cfg.rangeLimits(1) & allRange<=cfg.rangeLimits(2);
rangeAxis = allRange(keep);
rcAll = fft(adc.*w,cfg.nRangeFFT,1)/sum(w);
rcAll = rcAll.*exp(1i*2*pi*freq*meta.tCenter);
rc = rcAll(keep,:,:,:,:);
if cfg.calibrateChannels
    for m=1:cfg.nTx
        for r=1:cfg.nRx
            ch=(m-1)*cfg.nRx+r;
            rc(:,:,m,r,:) = rc(:,:,m,r,:)*exp(-1i*pi/180*cfg.channelPhaseDeg(ch));
        end
    end
end
% 此处保留 RVP；BP 中用 +j*pi*slope*tau^2 精确补偿。
end
