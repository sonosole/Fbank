module Fbank

using FFTW
export fbankparams
export initfbank

mutable struct fbankparams
    winLength::Int  # 分帧参数-帧长
    winShift::Int   # 分帧参数-帧移
    numBanks::Int   # 梅尔滤波参数-滤波器个数
    numFFT::Int     # 傅立叶变换点数
    alpha::Real     # 预加重系数，默认0.97
    fs::Int         # 采样率，一般16kHz
    maxfreq::Int    # 频域最大有效频率下标
    epsilon::Real   # 防止下溢系数
    function fbankparams(;
        winLength  = 256,
        winShift   = 128,
        numBanks   = 32,
        numFFT     = 256,
        alpha      = 0.97,
        fs         = 16000,
        epsilon    = 1e-6)
        maxfreq    = numFFT>>1
        new(winLength,winShift, numBanks,numFFT,alpha,fs, maxfreq,epsilon)
    end
end


function filterbanks(numBanks::Int, numFFT::Int, fs::Int)
    MAX   = numFFT>>1;                    # 正频率部分的下标最大值
    freq  = (0:(MAX-1))/MAX * fs/2;       # 下标映射到频率
    Fmel  = 2595*log10.(1 .+ freq/700);   # 频率映射到梅尔频率
    dFmel = Fmel[MAX]/(numBanks+1);       # 将Mel带平分
    bank  = zeros(numBanks, MAX);         # 滤波器的频域权重系数
    cFmel = 0.0;                          # 每个Mel频带的中心Mel频率
    for n = 1:numBanks
        cFmel = cFmel + dFmel
        for m = 1:MAX
            if ( Fmel[m] >= cFmel-dFmel ) && ( Fmel[m] <= cFmel+dFmel )
                bank[n,m] = 1.0 - abs( Fmel[m] - cFmel )/dFmel
            end
        end
    end
    return bank
end


function window(winLen::Int)
	hamming = 0.54 .- 0.46 .* cos.( 2*pi .* (0:(winLen-1))/(winLen-1) )
end


function filterwav(data, alpha)
    return (data[2:end] - alpha .* data[1:end-1])
end


function splitwav(data, win, winLength::Int, winShift::Int)
    numFrame = div(length(data)-winLength, winShift) + 1
    firstIds = (0:(numFrame-1)) .* winShift .+ 1  # 帧起始下标
    lasstIds = firstIds .+ (winLength - 1)        # 帧结束下标
    frames   = zeros(winLength, numFrame)
    for i = 1:numFrame
        frames[:,i] = data[firstIds[i]:lasstIds[i]] .* win
    end
    return frames, numFrame
end


function initfbank(params::fbankparams)
    winLength  = params.winLength
    winShift   = params.winShift
    numBanks   = params.numBanks
    numFFT     = params.numFFT
    maxfreq    = params.maxfreq
    alpha      = params.alpha
    fs         = params.fs
    epsilon    = params.epsilon

    winfunc    = window(winLength)
    melbank    = filterbanks(numBanks, numFFT, fs)

    function offlineFbank(wav)
        data  = filterwav(wav, alpha)                             # 滤波
        frames, n = splitwav(data, winfunc, winLength, winShift)  # 分帧

        tmp = fft(frames, 1)                    # 时域到频域,按列计算
        pow = abs2.(tmp[1:maxfreq,:])           # 功率谱,提取有用部分
        return log.(melbank * pow .+ epsilon)   # 对数梅尔功率谱
    end
    return offlineFbank
end


function Base.show(io::IO, f::fbankparams)
    println(io, "———————————————————————")
    println(io, "  winLength = $(f.winLength)")
    println(io, "   winShift = $(f.winShift)")
    println(io, "   numBanks = $(f.numBanks)")
    println(io, "     numFFT = $(f.numFFT)")
    println(io, "      alpha = $(f.alpha)")
    println(io, "         fs = $(f.fs)")
    println(io, "    maxfreq = $(f.maxfreq)")
    println(io, "    epsilon = $(f.epsilon)")
    println(io, "———————————————————————")
end


end
