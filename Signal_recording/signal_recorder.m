samples_per_chip = 5;
record_duration = 3.2; % seconds

GPS_L1_carrier_frequency = 1575.42 * 10^6;
chips_per_period = 1023;
periods_per_second = 1000;
samples_per_period = samples_per_chip * chips_per_period;
periods_per_bit = 20;
bits_per_second = 50;
sample_rate = samples_per_period * periods_per_second;
%% 
Adalm.Gain = 40; % dB
Adalm.Freq_Central = GPS_L1_carrier_frequency;
Adalm.SampleRate = sample_rate;
Adalm.FrameSize = samples_per_period * periods_per_bit;
Adalm.FramesNumber = round (record_duration * bits_per_second); 

burst_length = Adalm.FramesNumber * Adalm.FrameSize;
Model.TimeSim = burst_length / Adalm.SampleRate;
%% 
Model.Name = "Pluto_receiver.slx";
load_system (Model.Name);
Out_Model = sim (Model.Name);
samples = Out_Model.IQ (1 : burst_length); % for some reason number of samples in rec_IQ is one more than FramesNumber
samples = transpose (samples);
%% 
save ("recorded_samples.mat", "samples", "-v7.3");