function [list_sattelite, timing_synch, fineFreq] = acquisition (iqData, samplingFreq, varargin)

%% парметры принятого сигнала
% samplingFreq - частота дискретизации (5 MHz)

codeFreqBasis = 1.023e6; % Базовая частота C/A-кода
acqSearchBand = 6;      % Диапазон поиска по частоте ±5 kHz
acqSearchStep = 0.05; % Шаг поиска 
acqThreshold = 9.5;      % Порог обнаружения сигнала
M = 10; % Когерентный/некогерентный подход 

if nargin>1
    for k = 1:2:(length(varargin)-1)
        if varargin{k}=='T'
            acqThreshold = varargin{k+1};
        end
        if varargin{k}=='M'
            M = varargin{k+1}; 
        end
    end
 end
samplesPerCode = round(samplingFreq / (codeFreqBasis / 1023)); % всего отсчетов на C/A code

resample_CA_all = zeros(32, samplesPerCode);

for  PRN = 1:32
    ca_code = CAcode(PRN);
    
    t_resampled = (0:1/samplingFreq:(1023/codeFreqBasis - 1/samplingFreq));
    
    resampled_CA = zeros(1, length(t_resampled));
    for i = 1:length(t_resampled)
        chip_index = floor(t_resampled(i) * codeFreqBasis) + 1;
        if chip_index > 1023
            chip_index = 1023;
        end
        resampled_CA(i) = ca_code(chip_index);
    end
    resample_CA_all(PRN, :) = resampled_CA;
end

numberOfFrqBins =  2*acqSearchBand/acqSearchStep + 1;   

acq_corr = zeros(32, numberOfFrqBins, samplesPerCode);
frqBins = 1e3*(-acqSearchBand:acqSearchStep:acqSearchBand);
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi / samplingFreq;

list_sattelite = [];
freq_offset = []; 
timing_synch = [];
fineFreq = [];

for  PRN = 1:32
    for idx_M = 1:M
        for idx_fo = 1:length(frqBins)
            sigCarr = exp(1i*frqBins(idx_fo) * phasePoints);
            IQfreqDom1 = sigCarr .* iqData(1+(idx_M-1)*samplesPerCode:(idx_M)*samplesPerCode);
            acq_corr(PRN, idx_fo, :) =  squeeze(acq_corr(PRN, idx_fo, :)).' +  abs(cxcorr_fft(IQfreqDom1,  resample_CA_all(PRN, :)).^2);             
        end
    end
    acq_corr = abs(acq_corr);
    [peak(PRN),idx_peak] = (max(acq_corr(PRN, :, :),[],'all'));
    [row, col] = ind2sub(size(squeeze(acq_corr(1,:,:))), idx_peak); 
    if (peak(PRN)/mean(mean(acq_corr(PRN, :, :))) > acqThreshold) 
        list_sattelite = [list_sattelite PRN];
        freq_offset = [freq_offset -frqBins(row)];
        timing_synch = [timing_synch col];

        caCode = CAcode(PRN); % Генерация C/A-кода для PRN
        codeValueIndex = floor(( (1:10*samplesPerCode)) / (samplingFreq/codeFreqBasis));
        longCaCode = caCode((rem(codeValueIndex, 1023) + 1));

        % Удаление C/A-кода из сигнала
        startIdx = col; % Начальная фаза кода из acquisition
        endIdx = startIdx + 10*samplesPerCode - 1;
        if endIdx > length(iqData)
            endIdx = length(iqData);
            longCaCode = longCaCode(1:endIdx-startIdx+1);
        end
        xCarrier = iqData(startIdx:endIdx) .* longCaCode;
        
       % FFT 
        nfft = 4*2^nextpow2(length(xCarrier));
        fftResult = fft(xCarrier, nfft);
        [~, maxFftIdx] = max(abs(fftResult));
        
       polyn =  polyfit([1 2 3],[abs(fftResult(maxFftIdx-1)) abs(fftResult(maxFftIdx)) abs(fftResult(maxFftIdx+1))],2);
       top =  (-polyn(2)/(2*polyn(1))-2);
       maxFftIdx = maxFftIdx + top;

        % Frequency calculation
        if maxFftIdx > nfft/2
            fineFreq = [fineFreq (maxFftIdx - nfft - 1) * samplingFreq/nfft];
        else
            fineFreq = [fineFreq (maxFftIdx - 1) * samplingFreq/nfft];
        end
        
    end
end


end


function [ h ] = cxcorr_fft( a,b )
    if (length(a) < length(b))
        c = [ a zeros(1,length(b)-length(a)) ];
        d = b;
    else
        c = a;
        d = [ b zeros(1,length(a)-length(b)) ];
    end
    % calculate crosscorrelation
    e = fft(c);
    f = fft(d);
    g = e.*conj(f);
    h = (ifft(g));
end

function ca_code = CAcode(PRN)
% Генерация C/A-кода GPS для указанного PRN (от 1 до 37)
% Вход: PRN - номер спутника (1-37)
% Выход: ca_code - вектор C/A-кода (1023 элемента, значения ±1)


% Инициализация регистров сдвига (G1 и G2) единицами
g1 = ones(1, 10); 
g2 = ones(1, 10);   

% Таблица отводов для G2 (по спецификации IS-GPS-200)
g2_taps = [
    2, 6;    % PRN 1
    3, 7;    % PRN 2
    4, 8;    % PRN 3
    5, 9;    % PRN 4
    1, 9;    % PRN 5
    2, 10;   % PRN 6
    1, 8;    % PRN 7
    2, 9;    % PRN 8
    3, 10;   % PRN 9
    2, 3;    % PRN 10
    3, 4;    % PRN 11
    5, 6;    % PRN 12
    6, 7;    % PRN 13
    7, 8;    % PRN 14
    8, 9;    % PRN 15
    9, 10;   % PRN 16
    1, 4;    % PRN 17
    2, 5;    % PRN 18
    3, 6;    % PRN 19
    4, 7;    % PRN 20
    5, 8;    % PRN 21
    6, 9;    % PRN 22
    1, 3;    % PRN 23
    4, 6;    % PRN 24
    5, 7;    % PRN 25
    6, 8;    % PRN 26
    7, 9;    % PRN 27
    8, 10;   % PRN 28
    1, 6;    % PRN 29
    2, 7;    % PRN 30
    3, 8;    % PRN 31
    4, 9;    % PRN 32
];

tap1 = g2_taps(PRN, 1);
tap2 = g2_taps(PRN, 2);
ca_code = zeros(1, 1023); 

% Генерация 1023 чипов кода
for i = 1:1023
    
    g2_out = mod(g2(tap1) + g2(tap2), 2);
    
    % C/A-код = XOR(G1, G2)
    ca_code(i) = 1 - 2*mod(g1(10) + g2_out, 2); % Преобразуем в ±1
    
    % Сдвигаем регистры
    g1 = [mod(g1(10) + g1(3), 2), g1(1:9)];
    g2 = [mod(g2(10) + g2(9) + g2(8) + g2(6) + g2(3) + g2(2), 2), g2(1:9)];
end
ca_code = -ca_code;

end
