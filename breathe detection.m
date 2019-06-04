function m = fcn(Ei, sigma_1)
waveform = phased.RectangularWaveform('SampleRate',25e12,...
    'PulseWidth',5e-11,'OutputFormat','Pulses',...
    'NumPulses',1,'PRF',2e10);
target = phased.RadarTarget('Model','Nonfluctuating',...
    'MeanRCS',1,'OperatingFrequency',1e9);
targetpos = phased.Platform('InitialPosition',[0; 0; 1100],...
    'Velocity',[0; 0; 0]);%% distance between X-Ray tube and detector
antenna = phased.IsotropicAntennaElement(...
    'FrequencyRange',[5e6 5e10]);
transmitter = phased.Transmitter('PeakPower',5e3,'Gain',20,...
'InUseOutputPort',true);
transpos = phased.Platform('InitialPosition',[0;0;0],...
    'Velocity',[0;0;0]);
radiator = phased.Radiator('OperatingFrequency',1e9,'Sensor',antenna);
collector = phased.Collector('OperatingFrequency',1e9,'Sensor',antenna);
channel = phased.FreeSpace('SampleRate',waveform.SampleRate,...
    'OperatingFrequency',1e9,'TwoWayPropagation',false);
receiver = phased.ReceiverPreamp('Gain',0,'LossFactor',0,...
    'SampleRate',5e6,'NoiseFigure',5,...
    'EnableInputPort',true,'SeedSource','Property','Seed',1e3);
NumPulses = 100;
sig = waveform(); % get waveform 
transpos = transpos.InitialPosition; % get transmitter position
rxsig = zeros(length(sig),NumPulses);
% transmit and receive ten pulses
for n = 1:NumPulses
    % update target position
    [tgtpos,tgtvel] = targetpos(1/waveform.PRF);
    [tgtrng,tgtang] = rangeangle(tgtpos,transpos);
    %tpos(n) = tgtrng;
    [txsig,txstatus] = transmitter(sig); % transmit waveform
    txsig = radiator(txsig,tgtang); % radiate waveform toward target
    txsig = channel(txsig,transpos,tgtpos,[0;0;0],tgtvel); % propagate waveform to target
    txsig = target(txsig); % reflect the signal
    % propagate waveform from the target to the transmiter
    txsig = channel(txsig,tgtpos,transpos,tgtvel,[0;0;0]);
    txsig = collector(txsig,tgtang); % collect signal
    rxsig(:,n) = receiver(txsig,~txstatus); % receive the signal
end
prf = waveform.PRF;
fs = waveform.SampleRate;
fasttime = unigrid(0,1/fs,1/prf,'[)');
rangebins = (physconst('LightSpeed')*fasttime)/2;
probfa = 1e-10;
NoiseBandwidth = 5e6/2;
npower = noisepow(NoiseBandwidth,...
    receiver.NoiseFigure,receiver.ReferenceTemperature);
thresh = npwgnthresh(probfa,NumPulses,'noncoherent');
thresh = sqrt(npower*db2pow(thresh));
[pks,range_detect] = findpeaks(pulsint(rxsig,'noncoherent'),...
    'MinPeakHeight',thresh,'SortStr','descend');
range_estimate = rangebins(range_detect(1));
ts = rxsig(range_detect(1),:).';
[Pxx,F] = periodogram(ts,[],256,prf,'centered')
% plot(F,10*log10(Pxx))
% grid
% xlabel('Frequency ')
% ylabel('Power (dB)')
% title('Periodogram Spectrum Estimate')
[Y,I] = max(Pxx);
lambda = physconst('LightSpeed')/1e9;
tgtspeed = dop2speed(F(I)/2,lambda);
    
A= zeros(126,1)
for i = 1:1:126
    
    A(i)=F(i+128)
end
m = zeros(126,1)
impedance = zeros(126,1)
for i = 1:1:126

    w = 2*pi*A(i);
    E0 = 8.854187817e-12;
    a(i) = (Ei)^.2 + (sigma_1/w)^.2;
    b(i) = (sigma_1/w*E0)^.2;
    x(i) = 377/(a(i)^(0.25));
y(i) = (1 + (b(i)/a(i)))^(0.25);
mag(i) = x(i)/y(i)
m(i) = mag(i)
end
impedance(i) = m(i)
end
