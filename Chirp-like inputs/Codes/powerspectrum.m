function[Psd, freq, PsdManSmooth, freqbin] = powerspectrum(t,v)

% Computation of the power spectrum (Psd) for a signal v(t)
% [t] = msec
% [f] = Hz


RAD = 0;
            

vm = mean(v);
y = v-vm;

T = t(end);
L = length(y);
Nfft = 2^nextpow2(L);
Y = fft(y,Nfft)/(L/2);
Ys = abs(Y(1:Nfft/2));
%Ys = abs(Y(1:Nfft/2).^2);
% if you want to smooth the PSD uncomment the next two lines
% Ysmooth = smooth(Ys,'lowess',9);
% Ys = Ysmooth;
Fs = L/(T/1000);
freq = Fs/2*linspace(0,1,Nfft/2);
Psd = Ys;

freqbin = freq;
PsdManSmooth = Psd;

% Manual Smooth (PsdManSmooth) of the "raw" power spectrum Psd

fmax = 1000;
fgraphcut = find(freq>fmax,1);
binsize = 1; %in Hz
binsize = floor(binsize./(freq(2)-freq(1)));
if binsize<=0
    binsize=1;
end


freqbin = zeros(1,ceil(fgraphcut/binsize));
PsdManSmooth = zeros(1,ceil(fgraphcut/binsize));

for k=1:ceil(fgraphcut/binsize)
    freqbin(k) = freq(k*binsize);
    PsdManSmooth(k) = max(Ys(k*binsize:(k+1)*binsize-1));
end

if RAD == 1
    freqbin = freqbin*2*pi/1000; % from Hz to rad
end

% figure
% hold on
% %plot(freq(1:fgraphcut),Psd(1:fgraphcut),'o','linewidth',1);
% plot(freqbin,PsdManSmooth,'o');
% axis([0 fmax,-0.01,max(Ys)*1.2]);
% set(gca,'fontsize',20);
% xlabel('Freq.  [Hz]');
% ylabel('|Y(f)|');

% figure
% hold on
% %plot(freq(1:fgraphcut),Psd(1:fgraphcut),'o','linewidth',1);
% plot(freqbin,PsdManSmooth,'o');
% axis([0 fmax,-0.01,max(Ys)*1.2]);
% set(gca,'fontsize',20);
% xlabel('Freq.  [Hz]');
% ylabel('|Y(f)|')




