% 2019-08-30

clearvars;
% close all;

Fsin=@(t,f,Ain) Ain*sin(2*pi*f*t/1000-pi/2);
Fsyn=@(t,Ain,t0,taudec) Ain*exp(-(t-t0)/taudec);
Fsqw=@(t,f,Ain,DC) Ain*square(2*pi*f*t/1000-pi/2,DC);

C = 1;
gL = 0.25; %0.25; %0.1;%0.25; %0.25;
g1 =0.25; %5; %0.25;
tau1 =  100;% 1; %100

% g1=0; %for 1D

a = -gL/C;
b = -g1/C;
c = 1/tau1;
d = -1/tau1;

Ain = 1;
dt = 0.1; %0.01
Ncycles = 100; %5000 %100

% Generation of the ZAP-like regular and random funtions
rng(10)
RNDM = 0;
if RNDM == 0
    CycleFreq = 1:1:Ncycles;
elseif RNDM == 1
    CycleFreq = randperm(Ncycles,Ncycles);
end
% CycleFreq=CycleFreq(end:-1:1);


%%%%Selectively randomize
if(0)
    rand_times = 50;
    for i = 1:rand_times
        aux1= randi(50)+50;
        aux2= randi(50);
        repo1 = CycleFreq(aux1);
        CycleFreq(aux1) = CycleFreq(aux2);
        CycleFreq(aux2) = repo1;
    end
end

% CycleFreq = nonLinspace(1,100,Ncycles,'exp10');
% CycleFreq = nonLinspace(1,100,Ncycles,'cos');
% CycleFreq = nonLinspace(1,100,Ncycles,'log10');


CyclePer = 1000./CycleFreq;

t = 0:dt:CyclePer(1);
Izaplike = Fsin(t,CycleFreq(1),Ain);
tbdry(1) = 0;
jbdry(1) = 1;
tbdry(2) = t(end);
jbdry(2) = length(t);
for k=2:Ncycles
    tloc = dt:dt:CyclePer(k);   
    Floc = Fsin(tloc,CycleFreq(k),Ain);
    t = [t t(end)+tloc];
    Izaplike = [Izaplike Floc];
    tbdry(k+1) = t(end);
    jbdry(k+1) = length(t);
end
Tmax = t(end);

% Computation of the numerical solution

v = zeros(1,length(t));
w = zeros(1,length(t));

for j=1:length(t)-1
    k1v = -gL*v(j)-g1*w(j)+Izaplike(j);
    k1w = (v(j)-w(j))/tau1;
    av = v(j)+k1v*dt;
    aw = w(j)+k1w*dt;
    k2v = -gL*av-g1*aw+Izaplike(j+1);
    k2w = (av-aw)/tau1;
    v(j+1) = v(j)+(k1v+k2v)*dt/2;
    w(j+1) = w(j)+(k1w+k2w)*dt/2;
end

% Computation of the peaks and peak-times
% no trough and trough times computation (vmin-->vmax)

vmax = zeros(1,Ncycles);
jpeak = zeros(1,Ncycles);
tpeak = zeros(1,Ncycles);
Phase = zeros(1,Ncycles);
vmin = zeros(1,Ncycles);
jtrough = zeros(1,Ncycles);

wmax = zeros(1,Ncycles);
wmin = zeros(1,Ncycles);
for k=1:Ncycles
    [vmax(k),jpeak(k)] = max(v(jbdry(k):jbdry(k+1)));    
    tpeak(k) = tbdry(k)+jpeak(k)*dt;
    Phase(k) = (tpeak(k)-tbdry(k)-CyclePer(k)/2)*CycleFreq(k);
    [vmin(k),jtrough(k)] = min(v(jbdry(k):jbdry(k+1))); 
    [wmax(k),~] = max(w(jbdry(k):jbdry(k+1)));  
    [wmin(k),~] = min(w(jbdry(k):jbdry(k+1)));  
end
Phase = Phase*2*pi/1000;

[CycleFreqOrd,vmaxord,phaseord,vminord] = OrderingCycleFreq(CycleFreq,vmax,Phase,vmin);
% [Z,Phi,fres,Zmax,QZ] = Impedance2D(a,b,c,d,CycleFreqOrd);


% figure
% hold on
% plot(t,Izaplike,'-b','linewidth',2);
% axis([0 Tmax -1.2 1.2]);
% set(gca,'fontsize',24);
% xlabel('t');
% ylabel('');
% 
% 
% figure
% hold on
% plot(t,v,'-b','linewidth',2);
% plot(t,Izaplike-8,'-','Color',[.7 .7 .7],'linewidth',2);
% axis([0 Tmax -10 5]);
% set(gca,'fontsize',24);
% xlabel('t  [ms]');
% ylabel('');

%%%% FFT
Zfft = fft(v)./fft(Izaplike);
Phifft = -atan(imag(Zfft)./real(Zfft));
Zfft = abs(Zfft);
df = 1000/t(end);
f= 0:df:length(Zfft)*df-df;

% figure
% hold on
% plot(CycleFreqOrd,vmaxord,'ob','linewidth',2);
% plot(CycleFreqOrd,Z,'-b','linewidth',2);  
% plot(f,Zfft,'.g','linewidth',2);  
% axis([0 100 0 10]);
% set(gca,'fontsize',24);
% xlabel('f  [Hz]');
% ylabel('v_{max}  [mV]');
% 
% ff = CycleFreqOrd;
% figure
% hold on
% plot(CycleFreqOrd,Phi,'-b','linewidth',2);
% plot(CycleFreqOrd,Phase,'ob','linewidth',2);
% plot(f,Phifft,'.g','linewidth',2);  
% plot([0 max(ff)],[0 0],'--','Color',[.7 .7 .7]);
% plot([0 max(ff)],[pi/2 pi/2],'--','Color',[.7 .7 .7]);
% plot([0 max(ff)],[-pi/2 -pi/2],'--','Color',[.7 .7 .7]);
% plot([0 max(ff)],[2*pi 2*pi],'--','Color',[.7 .7 .7]);
% axis([0 ff(end) -pi/2-0.2 pi/2+0.2])
% set(gca,'fontsize',24);
% xlabel('f  [Hz]');
% ylabel('\Phi');
% title('Phase profile');
% xlim([0,100])

% Generation of the ZAP-like regular and random synaptic-like funtions

% Ain = 1*2.5;
Ain=1;
taurse_e = 1;     % AMPA
taudec_e = 5;     % AMPA
taurse_i = 0.2;   % GABA_A     
taudec_i = 10;    % GABA_A
Eex = 20;
Ein = -60;

taudec = taudec_e;

t = 0:dt:CyclePer(1);
Isyntrain = Fsyn(t,Ain,t(1),taudec);
for k=2:Ncycles
    tloc = dt:dt:CyclePer(k);  
    Floc = Fsyn(tloc,Ain,tloc(1),taudec);
    t = [t t(end)+tloc];
    Isyntrain = [Isyntrain Floc];
end

% Computation of the numerical solution

vs = zeros(1,length(t));
ws = zeros(1,length(t));

for j=1:length(t)-1
    k1v = -gL*vs(j)-g1*ws(j)+Isyntrain(j);
    k1w = (vs(j)-ws(j))/tau1;
    av = vs(j)+k1v*dt;
    aw = ws(j)+k1w*dt;
    k2v = -gL*av-g1*aw+Isyntrain(j+1);
    k2w = (av-aw)/tau1;
    vs(j+1) = vs(j)+(k1v+k2v)*dt/2;
    ws(j+1) = ws(j)+(k1w+k2w)*dt/2;
end


% Computation of the peaks and peak-times

vsmax = zeros(1,Ncycles);
jpeak = zeros(1,Ncycles);
tpeak = zeros(1,Ncycles);
Phases = zeros(1,Ncycles);
vsmin = zeros(1,Ncycles);
ttrough = zeros(1,Ncycles);
jtrough = zeros(1,Ncycles);
wsmax = zeros(1,Ncycles);
wsmin = zeros(1,Ncycles);
for k=1:Ncycles
    [vsmax(k),jpeak(k)] = max(vs(jbdry(k):jbdry(k+1)));    
    tpeak(k) = tbdry(k)+jpeak(k)*dt;
    Phases(k) = (tpeak(k)-tbdry(k)-CyclePer(k)/2)*CycleFreq(k);   
    [vsmin(k),jtrough(k)] = min(vs(jbdry(k):jbdry(k+1))); 
    [wsmax(k),~] = max(ws(jbdry(k):jbdry(k+1)));  
    [wsmin(k),~] = min(ws(jbdry(k):jbdry(k+1)));  
end
Phases = Phases*2*pi/1000;

[CycleFreqOrd,vsmaxord,phasesord,vsminord] = OrderingCycleFreq(CycleFreq,vsmax,Phases,vsmin);


%%%% FFT
Zffts = fft(vs)./fft(Isyntrain);
Phiffts = -atan(imag(Zffts)./real(Zffts));
Zffts = abs(Zffts);
df = 1000/t(end);
f= 0:df:length(Zffts)*df-df;

% figure
% hold on
% plot(t,Isyntrain,'-b','linewidth',2);
% axis([0 Tmax -1.2 1.2]);
% set(gca,'fontsize',24);
% xlabel('t');
% ylabel('');
% 
% figure
% hold on
% plot(t,vs,'-b','linewidth',2);
% plot(t,Isyntrain-8,'-','Color',[.7 .7 .7],'linewidth',2);
% axis([0 Tmax -10 5]);
% set(gca,'fontsize',24);
% xlabel('t  [ms]');
% ylabel('');
% 
% figure
% hold on
% plot(CycleFreqOrd,vsmaxord,'ob','linewidth',2);
% plot(CycleFreqOrd,vsminord,'ob','linewidth',2);
% plot(CycleFreqOrd,(vsmaxord-vsminord)/2,'or','linewidth',2);
% plot(f,Zffts,'.g','linewidth',2);  
% axis([0 100 -10 10]);
% set(gca,'fontsize',24);
% xlabel('f  [Hz]');
% ylabel('v_{max}  [mV]');
% 
% figure
% hold on
% % plot(CycleFreqOrd,Phi,'-b','linewidth',2);
% plot(CycleFreqOrd,Phases,'ob','linewidth',2);
% plot(f,Phiffts,'.g','linewidth',2);  
% plot([0 max(ff)],[0 0],'--','Color',[.7 .7 .7]);
% plot([0 max(ff)],[pi/2 pi/2],'--','Color',[.7 .7 .7]);
% plot([0 max(ff)],[-pi/2 -pi/2],'--','Color',[.7 .7 .7]);
% plot([0 max(ff)],[2*pi 2*pi],'--','Color',[.7 .7 .7]);
% axis([0 ff(end) -pi/2-0.2 pi/2+0.2])
% set(gca,'fontsize',24);
% xlabel('f  [Hz]');
% ylabel('\Phi');
% title('Phase profile');
% xlim([0,100])


% Generation of the ZAP-like square-wave funtions

DC = 50;
Ain = 1;

t = 0:dt:CyclePer(1);
Isqwtrain = Fsqw(t,CycleFreq(1),Ain,DC);
for k=2:Ncycles
    tloc = dt:dt:CyclePer(k);  
    Floc = Fsqw(tloc,CycleFreq(k),Ain,DC);   
    t = [t t(end)+tloc];
    Isqwtrain = [Isqwtrain Floc];
end

% Computation of the numerical solution

vsqw = zeros(1,length(t));
wsqw = zeros(1,length(t));

for j=1:length(t)-1
    k1v = -gL*vsqw(j)-g1*wsqw(j)+Isqwtrain(j);
    k1w = (vsqw(j)-wsqw(j))/tau1;
    av = vsqw(j)+k1v*dt;
    aw = wsqw(j)+k1w*dt;
    k2v = -gL*av-g1*aw+Isqwtrain(j+1);
    k2w = (av-aw)/tau1;
    vsqw(j+1) = vsqw(j)+(k1v+k2v)*dt/2;
    wsqw(j+1) = wsqw(j)+(k1w+k2w)*dt/2;
end

% Computation of the peaks and peak-times

vsqwmax = zeros(1,Ncycles);
jpeak = zeros(1,Ncycles);
tpeak = zeros(1,Ncycles);
Phasesqw = zeros(1,Ncycles);
vsqwmin = zeros(1,Ncycles);
ttrough = zeros(1,Ncycles);
jtrough = zeros(1,Ncycles);
wsqwmax = zeros(1,Ncycles);
wsqwmin = zeros(1,Ncycles);
for k=1:Ncycles
    [vsqwmax(k),jpeak(k)] = max(vsqw(jbdry(k):jbdry(k+1)));    
    tpeak(k) = tbdry(k)+jpeak(k)*dt;
    Phasesqw(k) = (tpeak(k)-tbdry(k)-CyclePer(k)/2)*CycleFreq(k);   
    [vsqwmin(k),jtrough(k)] = min(vsqw(jbdry(k):jbdry(k+1))); 
    [wsqwmax(k),~] = max(wsqw(jbdry(k):jbdry(k+1)));  
    [wsqwmin(k),~] = min(wsqw(jbdry(k):jbdry(k+1)));  
    
end
Phasesqw = Phasesqw*2*pi/1000;

[CycleFreqOrd,vsqwmaxord,phasesqword,vsqwminord] = OrderingCycleFreq(CycleFreq,vsqwmax,Phasesqw,vsqwmin);

%%%% FFT
Zfftsqw = fft(vsqw)./fft(Isqwtrain);
Phifftsqw = -atan(imag(Zfftsqw)./real(Zfftsqw));
Zfftsqw = abs(Zfftsqw);
df = 1000/t(end);
f= 0:df:length(Zfftsqw)*df-df;


% figure
% hold on
% plot(t,Isqwtrain,'-b','linewidth',2);
% axis([0 Tmax -1.2 1.2]);
% set(gca,'fontsize',24);
% xlabel('t');
% ylabel('');
% 
% figure
% hold on
% plot(t,vsqw,'-b','linewidth',2);
% plot(t,Isqwtrain-8,'-','Color',[.7 .7 .7],'linewidth',2);
% axis([0 Tmax -10 5]);
% set(gca,'fontsize',24);
% xlabel('t  [ms]');
% ylabel('');
% 
% figure
% hold on
% plot(CycleFreqOrd,vsqwmaxord,'ob','linewidth',2);
% plot(CycleFreqOrd,vsqwminord,'ob','linewidth',2);
% plot(CycleFreqOrd,(vsqwmaxord-vsqwminord)/2,'or','linewidth',2);
% plot(f,Zfftsqw,'.g','linewidth',2);  
% axis([0 100 -10 10]);
% set(gca,'fontsize',24);
% xlabel('f  [Hz]');
% ylabel('v_{max}  [mV]');
% 
% figure
% hold on
% % plot(CycleFreqOrd,Phi,'-b','linewidth',2);
% plot(CycleFreqOrd,Phasesqw,'ob','linewidth',2);
% plot(f,Phifftsqw,'.g','linewidth',2);  
% plot([0 max(ff)],[0 0],'--','Color',[.7 .7 .7]);
% plot([0 max(ff)],[pi/2 pi/2],'--','Color',[.7 .7 .7]);
% plot([0 max(ff)],[-pi/2 -pi/2],'--','Color',[.7 .7 .7]);
% plot([0 max(ff)],[2*pi 2*pi],'--','Color',[.7 .7 .7]);
% axis([0 ff(end) -pi/2-0.2 pi/2+0.2])
% set(gca,'fontsize',24);
% xlabel('f  [Hz]');
% ylabel('\Phi');
% title('Phase profile');
% xlim([0,100])

%  % Power spectrum density computation
%     
% [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,Izaplike);
% ZAP_Psd = PsdManSmooth;
% 
% [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,v);
% V_Psd = PsdManSmooth;
% 
% Z = V_Psd./ZAP_Psd;
% 
% figure
% hold on
% plot(freqbin,Z,'ob','linewidth',2);
% plot(freqbin,Z,'-b','linewidth',1);
% axis([0 200,-0.01,20]);
% set(gca,'fontsize',20);
% xlabel('Freq.  [Hz]');
% ylabel('Z');


figure; 
subplot(3,4,1)
plot(t,v,'-b','linewidth',2);
axis([0 Tmax -6 5]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('v [mV]');
title('voltage series')

subplot(3,4,2)
hold on
plot(CycleFreqOrd,vmaxord,'ob','linewidth',2);
plot(CycleFreqOrd,vminord,'ob','linewidth',2);
axis([0 100 -6 5]);
set(gca,'fontsize',24);
xlabel('f  [Hz]');
ylabel('v_{max}  [mV]');
title('envelopes')

subplot(3,4,3)
hold on;
plot(CycleFreqOrd,(vmaxord-vminord)/2,'ob','linewidth',2);
plot(f,Zfft,'.g','linewidth',2);  
axis([0 100 0 5]);
set(gca,'fontsize',24);
xlabel('f  [Hz]');
ylabel('Z');
legend('peak/trough','fft')
title('impedance')

subplot(3,4,4)
ff = CycleFreqOrd;
hold on
plot(CycleFreqOrd,Phase,'ob','linewidth',2);
plot(f,Phifft,'.g','linewidth',2);  
plot([0 max(ff)],[0 0],'--','Color',[.7 .7 .7]);
plot([0 max(ff)],[pi/2 pi/2],'--','Color',[.7 .7 .7]);
plot([0 max(ff)],[-pi/2 -pi/2],'--','Color',[.7 .7 .7]);
plot([0 max(ff)],[2*pi 2*pi],'--','Color',[.7 .7 .7]);
axis([0 ff(end) -pi/2-0.2 pi/2+0.2])
set(gca,'fontsize',24);
xlabel('f  [Hz]');
ylabel('\Phi');
title('phase');
xlim([0,100])

subplot(3,4,5)
plot(t,vs,'-b','linewidth',2);
axis([0 Tmax -6 5]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('v [mV]');

subplot(3,4,6)
hold on
plot(CycleFreqOrd,vsmaxord,'ob','linewidth',2);
plot(CycleFreqOrd,vsminord,'ob','linewidth',2);
axis([0 100 -6 5]);
set(gca,'fontsize',24);
xlabel('f  [Hz]');
ylabel('v_{max}  [mV]');

subplot(3,4,7)
hold on;
plot(CycleFreqOrd,(vsmaxord-vsminord)/2,'ob','linewidth',2);
plot(f,Zffts,'.g','linewidth',2);  
axis([0 100 0 5]);
set(gca,'fontsize',24);
xlabel('f  [Hz]');
ylabel('Z');

subplot(3,4,8)
hold on
plot(CycleFreqOrd,Phases,'ob','linewidth',2);
plot(f,Phiffts,'.g','linewidth',2);  
plot([0 max(ff)],[0 0],'--','Color',[.7 .7 .7]);
plot([0 max(ff)],[pi/2 pi/2],'--','Color',[.7 .7 .7]);
plot([0 max(ff)],[-pi/2 -pi/2],'--','Color',[.7 .7 .7]);
plot([0 max(ff)],[2*pi 2*pi],'--','Color',[.7 .7 .7]);
axis([0 ff(end) -pi/2-0.2 pi/2+0.2])
set(gca,'fontsize',24);
xlabel('f  [Hz]');
ylabel('\Phi');
xlim([0,100])

subplot(3,4,9)
plot(t,vsqw,'-b','linewidth',2);
axis([0 Tmax -6 5]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('v [mV]');

subplot(3,4,10)
hold on
plot(CycleFreqOrd,vsqwmaxord,'ob','linewidth',2);
plot(CycleFreqOrd,vsqwminord,'ob','linewidth',2);
axis([0 100 -6 5]);
set(gca,'fontsize',24);
xlabel('f  [Hz]');
ylabel('v_{max}  [mV]');

subplot(3,4,11)
hold on;
plot(CycleFreqOrd,(vsqwmaxord-vsqwminord)/2,'ob','linewidth',2);
plot(f,Zfftsqw,'.g','linewidth',2);  
axis([0 100 0 5]);
set(gca,'fontsize',24);
xlabel('f  [Hz]');
ylabel('Z');

subplot(3,4,12)
hold on;
plot(CycleFreqOrd,Phasesqw,'ob','linewidth',2);
plot(f,Phifftsqw,'.g','linewidth',2);  
plot([0 max(ff)],[0 0],'--','Color',[.7 .7 .7]);
plot([0 max(ff)],[pi/2 pi/2],'--','Color',[.7 .7 .7]);
plot([0 max(ff)],[-pi/2 -pi/2],'--','Color',[.7 .7 .7]);
plot([0 max(ff)],[2*pi 2*pi],'--','Color',[.7 .7 .7]);
axis([0 ff(end) -pi/2-0.2 pi/2+0.2])
set(gca,'fontsize',24);
xlabel('f  [Hz]');
ylabel('\Phi');
xlim([0,100])


% smoothedSignal = conv(vsqw, ones(1,2), 'same'); % or whatever
% newSignal = vsqw; % Initialize
% % Replace region of interest with smoothed values.
% newSignal(1:end) = smoothedSignal(1:end); % You determine the indexes in some way.
% plot(vsqw)
% hold on;
% plot(newSignal)

% results = [vmax; wmax; vmin; wmin]';
% save('/Users/rodrigo/Google Drive/project_pos_doc/codigos_impedance/Resonance_Cycles_randomize/linear_high_freq/envelopes_data/freqPlane_zap_linear_shuffled.dat','results','-ascii')
% results = [vsmax; wsmax; vsmin; wsmin]';
% save('/Users/rodrigo/Google Drive/project_pos_doc/codigos_impedance/Resonance_Cycles_randomize/linear_high_freq/envelopes_data/freqPlane_s_linear_shuffled.dat','results','-ascii')
% results = [vsqwmax; wsqwmax; vsqwmin; wsqwmin]';
% save('/Users/rodrigo/Google Drive/project_pos_doc/codigos_impedance/Resonance_Cycles_randomize/linear_high_freq/envelopes_data/freqPlane_sqw_linear_shuffled.dat','results','-ascii')

if(0) %no shuffled
    f_zap_sh = CycleFreqOrd;
    z_zap_sh = (vmaxord-vminord)/2;
    f_zap_sh_fft = f;
    z_zap_sh_fft = Zfft;

    results = [f_zap_sh'  z_zap_sh'];
    save('zap_pkth.dat','results','-ascii')

    results = [f_zap_sh_fft'  z_zap_sh_fft'];
    save('zap_fft.dat','results','-ascii')

    f_exp_sh = CycleFreqOrd;
    z_exp_sh = (vsmaxord-vsminord)/2;
    f_exp_sh_fft = f;
    z_exp_sh_fft = Zffts;

    results = [f_exp_sh'  z_exp_sh'];
    save('exp_pkth.dat','results','-ascii')

    results = [f_exp_sh_fft'  z_exp_sh_fft'];
    save('exp_fft.dat','results','-ascii')

    f_sqw_sh = CycleFreqOrd;
    z_sqw_sh = (vsqwmaxord-vsqwminord)/2;
    f_sqw_sh_fft = f;
    z_sqw_sh_fft = Zfftsqw;

    results = [f_sqw_sh'  z_sqw_sh'];
    save('sqw_pkth.dat','results','-ascii')

    results = [f_sqw_sh_fft'  z_sqw_sh_fft'];
    save('sqw_fft.dat','results','-ascii')
end

if(0) %no shuffled   
    save('v.dat','v','-ascii');
    save('vs.dat','vs','-ascii');
    save('vsqw.dat','vsqw','-ascii');
    
    save('v_sh.dat','v','-ascii');
    save('vs_sh.dat','vs','-ascii');
    save('vsqw_sh.dat','vsqw','-ascii');
    
    save('vmax.dat','vmaxord','-ascii');
    save('vsmax.dat','vsmaxord','-ascii');
    save('vsqwmax.dat','vsqwmaxord','-ascii');
    save('vmin.dat','vminord','-ascii');
    save('vsmin.dat','vsminord','-ascii');
    save('vsqwmin.dat','vsqwminord','-ascii');
    
    save('vmax_sh.dat','vmaxord','-ascii');
    save('vsmax_sh.dat','vsmaxord','-ascii');
    save('vsqwmax_sh.dat','vsqwmaxord','-ascii');
    save('vmin_sh.dat','vminord','-ascii');
    save('vsmin_sh.dat','vsminord','-ascii');
    save('vsqwmin_sh.dat','vsqwminord','-ascii');
    
end

if(0) %shuffled
    f_zap_sh = CycleFreqOrd;
    z_zap_sh = (vmaxord-vminord)/2;
    f_zap_sh_fft = f;
    z_zap_sh_fft = Zfft;

    results = [f_zap_sh'  z_zap_sh'];
    save('zap_sh_pkth.dat','results','-ascii')

    results = [f_zap_sh_fft'  z_zap_sh_fft'];
    save('zap_sh_fft.dat','results','-ascii')

    f_exp_sh = CycleFreqOrd;
    z_exp_sh = (vsmaxord-vsminord)/2;
    f_exp_sh_fft = f;
    z_exp_sh_fft = Zffts;

    results = [f_exp_sh'  z_exp_sh'];
    save('exp_sh_pkth.dat','results','-ascii')

    results = [f_exp_sh_fft'  z_exp_sh_fft'];
    save('exp_sh_fft.dat','results','-ascii')

    f_sqw_sh = CycleFreqOrd;
    z_sqw_sh = (vsqwmaxord-vsqwminord)/2;
    f_sqw_sh_fft = f;
    z_sqw_sh_fft = Zfftsqw;

    results = [f_sqw_sh'  z_sqw_sh'];
    save('sqw_sh_pkth.dat','results','-ascii')

    results = [f_sqw_sh_fft'  z_sqw_sh_fft'];
    save('sqw_sh_fft.dat','results','-ascii')
end

% szap = (fft(Izaplike));
% ssyn = (fft(Isyntrain));
% ssqw = (fft(Isqwtrain));
% 
% svzap = (fft(v));
% svsyn = (fft(vs));
% svsqw = (fft(vsqw));
% figure; 
% subplot(3,2,1)
% pps = -atan(imag(szap)./real(szap));
% ppv = -atan(imag(svzap)./real(svzap));
% plot(f,pps); hold on; 
% plot(f,ppv); hold on; 
% plot(f,real(svzap));hold on; 
% plot(f,real(szap)); hold on; 
% % plot(f,imag(svzap));
% xlim([0,100])
% subplot(3,2,3)
% semilogy(f,ssyn); hold on; 
% semilogy(f,svsyn);
% xlim([0,100])
% subplot(3,2,5)
% semilogy(f,ssqw); hold on;
% semilogy(f,svsqw);
% xlim([0,100])
% 
% subplot(3,2,2)
% plot(f,(svzap)./(szap),'.g','linewidth',2);  axis([0 100 0 5]);
% subplot(3,2,4)
% plot(f,Zffts,'.g','linewidth',2);  axis([0 100 0 5]);
% subplot(3,2,6)
% plot(f,Zfftsqw,'.g','linewidth',2);  axis([0 100 0 5]);
