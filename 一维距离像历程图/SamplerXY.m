function interpolatedValue = SamplerXY(x, y, freqStart, fRange, freqSamples, angStart, angRange, angSamples, Sampler)
    % 计算极坐标
    freq = sqrt(x^2 + y^2);
    ang = rad2deg(atan2(y, x));
    
    % 归一化采样索引
    freqSample = clamp((freq - freqStart) / fRange * (freqSamples - 1), 0, freqSamples - 1) + 1;
    angSample  = clamp((ang - angStart) / angRange * (angSamples - 1), 0, angSamples - 1) + 1;
    
    % 确保索引是整数
    freqSampleL = floor(freqSample);
    freqSampleU = ceil(freqSample);
    angSampleL  = floor(angSample);
    angSampleU  = ceil(angSample);
    
    % 计算在区间内的相对位置
    freqSampleD = freqSample - freqSampleL;
    angSampleD  = angSample - angSampleL;
 
    % 双线性插值
    ll = freqSampleD * angSampleD;
    ul = (1 - freqSampleD) * angSampleD;
    uu = (1 - freqSampleD) * (1 - angSampleD);
    lu = freqSampleD * (1 - angSampleD);
 
    % 返回：四个邻近采样点的加权平均值
    interpolatedValue = uu * Sampler(freqSampleL, angSampleL) + ...
                        ul * Sampler(freqSampleL, angSampleU) + ...
                        lu * Sampler(freqSampleU, angSampleL) + ...
                        ll * Sampler(freqSampleU, angSampleU);
end