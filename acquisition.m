clc; 
close all; 
clear all;

%% парметры принятого сигнала
samplingFreq = 5e6;      % Частота дискретизации 5 MHz
codeFreqBasis = 1.023e6; % Базовая частота C/A-кода
acqSearchBand = 6;      % Диапазон поиска по частоте ±5 kHz
acqSearchStep = 0.05; % Шаг поиска 
acqThreshold = 4;      % Порог обнаружения сигнала
M = 3; % Когерентный/некогерентный подход 
%% считаем 1 мс с файла 
samplesPerCode = round(samplingFreq / (codeFreqBasis / 1023)); % всего отсчетов на C/A code
samples100ms = 100 * samplesPerCode; % Отсчетов в 100 мс (С/A code имеет длительность в 1 мс)

filename = 'GPS1min.bin';
fid = fopen(filename, 'rb');
if fid == -1
    error('Не удалось открыть файл %s', filename);
end

% Пропускаем первые 100 мс (умножаем на 2, так как I и Q чередуются)
fseek(fid, samples100ms * 2, 'bof');

% Читаем следующие 1*M мс для обработки (I и Q чередуются)
samplesToRead = M * samplesPerCode * 2;
data = fread(fid, samplesToRead, 'int8').';
fclose(fid);

% Преобразование в комплексный сигнал (I + jQ)
if mod(length(data), 2) ~= 0
    data = data(1:end-1); % Удаляем последний отсчет, если нечетное количество
end
iqData = data(1:2:end) + 1i * data(2:2:end);

iqData = iqData / 128; % нормализация 

%%  C/A code table 
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

%% acquizition 

numberOfFrqBins =  2*acqSearchBand/acqSearchStep + 1;   

acq_corr = zeros(32, numberOfFrqBins, samplesPerCode);
frqBins = 1e3*(-acqSearchBand:acqSearchStep:acqSearchBand);
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi / samplingFreq;

list_sattelite = [];
freq_offset = []; 
timing_synch = [];
peak_list = [];
fineFreq = [];

for  PRN = 1:32
    for idx_M = 1:M
        for idx_fo = 1:length(frqBins)
            sigCarr = exp(1i*frqBins(idx_fo) * phasePoints);
            IQfreqDom1 = sigCarr .* iqData(1+(idx_M-1)*samplesPerCode:(idx_M)*samplesPerCode);
            acq_corr(PRN, idx_fo, :) =  squeeze(acq_corr(PRN, idx_fo, :)).' +  abs(cxcorr_fft(IQfreqDom1,  resample_CA_all(PRN, :)));
        end
    end
    [peak(PRN),idx_peak] = (max(acq_corr(PRN, :, :),[],'all'));
    peak_list = [peak_list peak(PRN)/mean(mean(acq_corr(PRN, :, :)))];
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

 real_fo = [2.527236938476563e+03 5.245208740234375e+02 -2.098083496093750e+03 -2.822875976562500e+03 5.626678466796875e+02 2.269744873046875e+03];
error = (fineFreq - real_fo)./real_fo


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
