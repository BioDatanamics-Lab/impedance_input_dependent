% 2021-02-04

clearvars;
close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightgray = [.7 .7 .7];
darkgray = [.3 .3 .3];
darkgray2 = [.1 .1 .1];

Fsyn=@(t,Ain,t0,taudec) Ain*exp(-(t-t0)/taudec);

SPKSE = 3;
        % 1: Uniform
        % 2: Poisson
        % 3: Poisson (ISI>ISI_min)

C = 1;
gL = 0.25;
g1 = 0.25;
tau1 = 100;

% gL = 0.05;
% g1 = 0.3;

% g1 = 0;
%

C = 1;
gL = 0.1;
g1 = 0;
tau1 = 100;

D = 1;

taurse_e = 1;     % AMPA
taudec_e = 5;     % AMPA
taurse_i = 0.2;   % GABA_A     
taudec_i = 10;    % GABA_A
Eex = 20;
Ein = -60;

taudec_e = 25;

% Eigenvalues of the autonomous system 

a = -gL/C;
b = -g1/C;
c = 1/tau1;
d = -1/tau1;

[r,mu,fnat] = Eigenvalues2D(a,b,c,d);

% Impedance profile

freq = 1:1:100;
[Zanl,Phi,fres,Zmax,QZ] = Impedance2D(a,b,c,d,freq);

% Time definitions

Tmax = 1000000;
dt = 0.1;
t = 0:dt:Tmax;

Tmaxx = Tmax;

% White (Gaussian) noise

eta = randn(1,length(t));

% Spike-time: generation 
    
if SPKSE == 1
    Freq = 40;              
    [tspk,ISIpre] = Uniform(Freq,Tmax,dt);
elseif SPKSE == 2
    Freq = 1000;
    [tspk,ISIpre] = Poisson(Freq,Tmax,dt);            
elseif SPKSE == 3
    Freq = 1000;
    ISImin = 0.2;
    [tspk,ISIpre,cnt] = PoissonMin(Freq,Tmax,dt,ISImin);
end
jspk = floor(tspk/dt);  
tspk = t(jspk);

%t = 0:dt:Tmax-cnt;

binsize = 10;
Edges = 0:binsize:4000;
[Hpre,Edgespre]=histcounts(ISIpre,Edges);
Hmaxpre = max(Hpre);
%Hmax = max(Hmaxpre,Hmaxpost);

figure
hold on
histogram(ISIpre,Edges);
%histogram(ISIpre,Edges,'Facecolor','b');
axis([0 Edges(end) 0 Hmaxpre+10]);
set(gca,'fontsize',24);
xlabel('ISI (pre)')

% figure
% hold on
% plot(tspk,'ob');
% axis([0 length(tspk) 0 tspk(end)]);
% set(gca,'fontsize',24);
% xlabel('Spike #')
% ylabel('t_{spk}');

Ain = 1;
Isyntrain = zeros(1,length(t));
for k=1:length(tspk)-1
    Isyntrain(jspk(k):jspk(k+1)-1) = Fsyn(t(jspk(k):jspk(k+1)-1),Ain,tspk(k),taudec_e);
end
    
% figure
% hold on
% plot(t,Isyntrain,'-b','linewidth',2);
% axis([0 Tmax 0 1.2]);
% set(gca,'fontsize',24);
% xlabel('t')
% ylabel('I_{syn}');


% Excitatory synaptic-like current-based inputs

Gsyn = 1;

% Computation of the numerical solution

vs = zeros(1,length(t));
ws = zeros(1,length(t));

for j=1:length(t)-1
    k1v = -gL*vs(j)-g1*ws(j)+Gsyn*Isyntrain(j);
    k1w = (vs(j)-ws(j))/tau1;
    av = vs(j)+k1v*dt;
    %av = av+sqrt(2*D*dt)*eta(j);
    aw = ws(j)+k1w*dt;
    k2v = -gL*av-g1*aw+Gsyn*Isyntrain(j+1);
    k2w = (av-aw)/tau1;
    vs(j+1) = vs(j)+(k1v+k2v)*dt/2;
    %vs(j+1) = vs(j+1)+sqrt(2*D*dt)*eta(j);
    ws(j+1) = ws(j)+(k1w+k2w)*dt/2;
end

vslin = vs;
wslin = ws;

% figure
% hold on
% plot(t,vs,'-b','linewidth',2);
% plot(t,Isyntrain-3,'-','Color',[.7 .7 .7],'linewidth',2);
% %plot(tsepeak,vsemax,'or')
% %plot(tsetrough,vsemin,'or')
% axis([0 Tmax -4 5]);
% set(gca,'fontsize',24);
% xlabel('t  [ms]');
% ylabel('V');
% legend('Current-based');


% % Computation of the peaks and peak-times
% 
% [tpeakin,maxin,contpeakin] = PeaksOsc(Isyntrain,t,dt,Tmax,dt);
% [tspeak,vsmax,cntpeak] = PeaksOsc(vs,t,dt,Tmax,dt);
% [tstrough,vsmin,cnttrough] = TroughsOsc(vs,t,dt,Tmax,dt);
% CycleFreqOutputPeak = 1000./diff(tstrough);
% [CycleFreqOrdPeak,vsmaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsmax(2:end));
% CycleFreqOutputTrough = 1000./diff(tspeak);
% [CycleFreqOrdTrough,vsminord] = OrderingCycleFreqTrough(CycleFreqOutputTrough,vsmin(1:end-1)); 
% 
% figure
% hold on
% plot(CycleFreqOrdPeak,vsmaxord,'ob','linewidth',2);
% plot(CycleFreqOrdTrough,vsminord,'ob','linewidth',2);
% axis([0 100 -2 4]);
% set(gca,'fontsize',24);
% xlabel('f  [Hz]');
% ylabel('v_{max}  [mV]');



%10: 150
%100: 90

% Excitatory synaptic-like conductance-based inputs

Gsyncb = 1;
    
% Excitatory 

Esyn = 1;

% Computation of the numerical solution

vscb = zeros(1,length(t));
wscb = zeros(1,length(t));

for j=1:length(t)-1
    k1v = -gL*vscb(j)-g1*wscb(j)-Gsyncb*Isyntrain(j)*(vscb(j)-Esyn);
    k1w = (vscb(j)-wscb(j))/tau1;
    av = vscb(j)+k1v*dt;
    av = av+sqrt(2*D*dt)*eta(j);
    aw = wscb(j)+k1w*dt;
    k2v = -gL*av-g1*aw-Gsyncb*Isyntrain(j+1)*(av-Esyn);
    k2w = (av-aw)/tau1;
    vscb(j+1) = vscb(j)+(k1v+k2v)*dt/2;
    vscb(j+1) = vscb(j+1)+sqrt(2*D*dt)*eta(j);
    wscb(j+1) = wscb(j)+(k1w+k2w)*dt/2;
end

vscbe = vscb;
wscbe = wscb;

vscbelin = vscbe;
wscbelin = wscbe;

% % Computation of the peaks and peak-times
%     
% vscbmax = zeros(1,Ncycles);
% jpeak = zeros(1,Ncycles);
% tpeak = zeros(1,Ncycles);
% Phasescb = zeros(1,Ncycles);
% vscbmin = zeros(1,Ncycles);
% ttrough = zeros(1,Ncycles);
% jtrough = zeros(1,Ncycles);
% for k=1:Ncycles
%     [vscbmax(k),jpeak(k)] = max(vscb(jbdry(k):jbdry(k+1)));    
%     tpeak(k) = tbdry(k)+jpeak(k)*dt;
%     Phasescb(k) = (tpeak(k)-tbdry(k)-CyclePer(k)/2)*CycleFreq(k);   
%     [vscbmin(k),jtrough(k)] = min(vscb(jbdry(k):jbdry(k+1))); 
%     ttrough(k) = tbdry(k)+jtrough(k)*dt;
% end
% Phasescb = Phasescb*2*pi/1000;

% Izaplikecbe = -(Gsyncb*Isyntrain.*(vscbe-Esyn));
% figure
% hold on
% plot(t,vscbe,'-b','linewidth',2);
% %plot(t,-Isyntrain.*(vscb-Eex)/max(-Isyntrain.*(vscb-Eex))-2,'-','Color',[.7 .7 .7],'linewidth',2);
% %plot(tscbepeak,vscbmax,'or','linewidth',2)
% %plot(tscbetrough,vscbmin,'or','linewidth',2)
% %plot(t,Izaplikecbe/max(Izaplikecbe)-2,'-','Color',[.7 .7 .7],'linewidth',2);
% plot(t,Izaplikecbe-2,'-','Color',[.7 .7 .7],'linewidth',2);
% axis([0 Tmax -3 5]);
% set(gca,'fontsize',24);
% xlabel('t  [ms]');
% ylabel('');
% legend('Conductance-based');

% White noise computations

% White (Gaussian) noise

eta = randn(1,length(t));

% Computation of the numerical solution

v = zeros(1,length(t));
w = zeros(1,length(t));
v(1) = 0;
w(1) = 0;

for j=1:length(t)-1
    k1v = -gL*v(j)-g1*w(j);
    k1w = (v(j)-w(j))/tau1;
    av = v(j)+k1v*dt;
    av = av+sqrt(2*D*dt)*eta(j);
    aw = w(j)+k1w*dt;
    k2v = -gL*av-g1*aw;
    k2w = (av-aw)/tau1;
    v(j+1) = v(j)+(k1v+k2v)*dt/2;
    v(j+1) = v(j+1)+sqrt(2*D*dt)*eta(j);
    w(j+1) = w(j)+(k1w+k2w)*dt/2;
end

vn = v;
wn = w;

vnlin = vn;
wnlin = wn;

% figure
% hold on
% plot(t,vn,'-b','linewidth',2);
% plot(t,eta-2,'-','Color',[.7 .7 .7],'linewidth',2);
% axis([0 Tmax -3 5]);
% set(gca,'fontsize',24);
% xlabel('t  [ms]');
% ylabel('');
% legend('White noise');


% Power spectrum density computation

[PSD,Freq,PsdManSmooth,freqbin] = powerspectrum(t,Isyntrain);
ZAP_Psd = PsdManSmooth;

[PSD,Freq, PsdManSmooth,freqbin] = powerspectrum(t,vslin);
Vs_Psd = PsdManSmooth;

Zs = Vs_Psd./ZAP_Psd;


[PSD,Freq, PsdManSmooth,freqbin] = powerspectrum(t,vscbelin);
Vscb_Psd = PsdManSmooth;

Zscb = Vscb_Psd./ZAP_Psd;

[PSD,Freq,PsdManSmooth,freqbin] = powerspectrum(t,eta);
ZAPn_Psd = PsdManSmooth;

[PSD,Freq, PsdManSmooth,freqbin] = powerspectrum(t,vnlin);
Vn_Psd = PsdManSmooth;


Zn = Vn_Psd./ZAPn_Psd;


figure
hold on
plot(ZAP_Psd,'ob','linewidth',1)
set(gca,'fontsize',20);
xlabel('Freq.  [Hz]');
ylabel('');

figure
hold on
plot(-100,-100,'ob','linewidth',2);
plot(-100,-100,'or','linewidth',2);
plot(-100,-100,'og','linewidth',2);
plot(freqbin,Zs/Gsyn,'ob','linewidth',1);
plot(freq,Zanl,'-r','linewidth',2);
plot(freqbin,Zn,'-g','linewidth',2);
axis([0 100,-0.01,5]);
set(gca,'fontsize',20);
xlabel('Freq.  [Hz]');
ylabel('Z');
legend('Z_{Poisson}','Z_{ANL}');



figure
hold on
plot(-100,-100,'ob','linewidth',2);
plot(-100,-100,'or','linewidth',2);
plot(-100,-100,'og','linewidth',2);
plot(freqbin,Vs_Psd,'ob','linewidth',1);
plot(freqbin,smooth(Vs_Psd,'moving',13),'-b','linewidth',2);
plot(freqbin,Vscb_Psd,'or','linewidth',1);
plot(freqbin,smooth(Vscb_Psd,'moving',13),'-r','linewidth',2);
plot(freqbin,Vn_Psd,'og','linewidth',1);
plot(freqbin,smooth(Vn_Psd,'moving',13),'-g','linewidth',2);
axis([0 100,-0.01,0.1]);
set(gca,'fontsize',24);
xlabel('Freq.  [Hz]');
ylabel('V_{PSD}');
legend('Current-based','Conductance-based','White noise');


figure
hold on
plot(-100,-100,'ob','linewidth',2);
plot(-100,-100,'or','linewidth',2);
plot(-100,-100,'og','linewidth',2);
plot(freqbin,Vs_Psd,'ob','linewidth',1);
plot(freqbin,smooth(Vs_Psd,'moving',13),'-b','linewidth',2);
plot(freqbin,Vscb_Psd,'or','linewidth',1);
plot(freqbin,smooth(Vscb_Psd,'moving',13),'-r','linewidth',2);
plot(freqbin,Vn_Psd,'og','linewidth',1);
plot(freqbin,smooth(Vn_Psd,'moving',13),'-g','linewidth',2);
plot(freqbin,smooth(Vs_Psd,'moving',13)*Vn_Psd(1)/Vs_Psd(1),'--b','linewidth',2);
plot(freqbin,smooth(Vscb_Psd,'moving',13)*Vn_Psd(1)/Vscb_Psd(1),'--r','linewidth',2);
%plot(freqbin,smooth(Vs_Psd,'moving',13)*Vn_Psd(1)/Vs_Psd(1),'-','Color',lightblueish,'linewidth',2);
%plot(freqbin,smooth(Vscb_Psd,'moving',13)*Vn_Psd(1)/Vscb_Psd(1),'-','Color',lightcoral,'linewidth',2);
axis([0 100,0,0.1]);
set(gca,'fontsize',24);
xlabel('Freq.  [Hz]');
ylabel('V_{PSD}');
legend('Current-based','Conductance-based','White noise');

ExIn =  0;
if ExIn == 1
    
    Gsyni = Gsyn;
    Gsyncbi = Gsyncb; 
    Esyni = -0.5;

    % Spike-time: generation 

    if SPKSE == 1
        Freq = 40;              
        [tspk,ISIpre] = Uniform(Freq,Tmax,dt);
    elseif SPKSE == 2
        Freq = 100;
        [tspk,ISIpre] = Poisson(Freq,Tmax,dt);            
    elseif SPKSE == 3
        Freq = 500;
        ISImin = 0.2;
        [tspk,ISIpre,cnt] = PoissonMin(Freq,Tmax,dt,ISImin);
    end
    jspk = floor(tspk/dt);  
    tspk = t(jspk);
    
    Ain = 1;
    Isyntraini = zeros(1,length(t));
    for k=1:length(tspk)-1
        Isyntraini(jspk(k):jspk(k+1)-1) = Fsyn(t(jspk(k):jspk(k+1)-1),Ain,tspk(k),taudec_i);
    end
    
    % Inhibitory synaptic-like current-based inputs

    % Computation of the numerical solution

    vs = zeros(1,length(t));
    ws = zeros(1,length(t));

    for j=1:length(t)-1
        k1v = -gL*vs(j)-g1*ws(j)+Gsyn*Isyntrain(j)-Gsyni*Isyntraini(j);
        k1w = (vs(j)-ws(j))/tau1;
        av = vs(j)+k1v*dt;
        %av = av+sqrt(2*D*dt)*eta(j);
        aw = ws(j)+k1w*dt;
        k2v = -gL*av-g1*aw+Gsyn*Isyntrain(j+1)-Gsyni*Isyntraini(j+1);
        k2w = (av-aw)/tau1;
        vs(j+1) = vs(j)+(k1v+k2v)*dt/2;
        %vs(j+1) = vs(j+1)+sqrt(2*D*dt)*eta(j);
        ws(j+1) = ws(j)+(k1w+k2w)*dt/2;
    end

    vslinei = vs;
    wslinei = ws;

    % Inhibitory synaptic-like conductance-based inputs

    % Computation of the numerical solution

    vscb = zeros(1,length(t));
    wscb = zeros(1,length(t));

    for j=1:length(t)-1
        k1v = -gL*vscb(j)-g1*wscb(j)-Gsyncb*Isyntrain(j)*(vscb(j)-Esyn)-Gsyncbi*Isyntraini(j)*(vscb(j)-Esyni);
        k1w = (vscb(j)-wscb(j))/tau1;
        av = vscb(j)+k1v*dt;
        av = av+sqrt(2*D*dt)*eta(j);
        aw = wscb(j)+k1w*dt;
        k2v = -gL*av-g1*aw-Gsyncb*Isyntrain(j+1)*(av-Esyn)-Gsyncbi*Isyntraini(j+1)*(av-Esyni);
        k2w = (av-aw)/tau1;
        vscb(j+1) = vscb(j)+(k1v+k2v)*dt/2;
        vscb(j+1) = vscb(j+1)+sqrt(2*D*dt)*eta(j);
        wscb(j+1) = wscb(j)+(k1w+k2w)*dt/2;
    end

    vscbei = vscb;
    wscbei = wscb;
    
    vscbeilin = vscbei;
    wscbeilin = wscbei;
    
    % Power spectrum density computation

    [PSD,Freq, PsdManSmooth,freqbin] = powerspectrum(t,vslinei);
    Vsei_Psd = PsdManSmooth;

    [PSD,Freq, PsdManSmooth,freqbin] = powerspectrum(t,vscbeilin);
    Vscbei_Psd = PsdManSmooth;

    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(-100,-100,'og','linewidth',2);
    plot(freqbin,Vsei_Psd,'ob','linewidth',1);
    plot(freqbin,smooth(Vsei_Psd,'moving',13),'-b','linewidth',2);
    plot(freqbin,Vscbei_Psd,'or','linewidth',1);
    plot(freqbin,smooth(Vscbei_Psd,'moving',13),'-r','linewidth',2);
    plot(freqbin,Vn_Psd,'og','linewidth',1);
    plot(freqbin,smooth(Vn_Psd,'moving',13),'-g','linewidth',2);
    axis([0 100,-0.01,0.1]);
    set(gca,'fontsize',24);
    xlabel('Freq.  [Hz]');
    ylabel('V_{PSD}');
    legend('Current-based','Conductance-based','White noise');
    
    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(-100,-100,'og','linewidth',2);
    plot(freqbin,Vsei_Psd,'ob','linewidth',1);
    plot(freqbin,smooth(Vsei_Psd,'moving',13),'-b','linewidth',2);
    plot(freqbin,Vscbei_Psd,'or','linewidth',1);
    plot(freqbin,smooth(Vscbei_Psd,'moving',13),'-r','linewidth',2);
    plot(freqbin,Vn_Psd,'og','linewidth',1);
    plot(freqbin,smooth(Vn_Psd,'moving',13),'-g','linewidth',2);
    plot(freqbin,smooth(Vsei_Psd,'moving',13)*Vn_Psd(1)/Vsei_Psd(1),'--b','linewidth',2);
    plot(freqbin,smooth(Vscbei_Psd,'moving',13)*Vn_Psd(1)/Vscbei_Psd(1),'--r','linewidth',2);
    %plot(freqbin,smooth(Vs_Psd,'moving',13)*Vn_Psd(1)/Vs_Psd(1),'-','Color',lightblueish,'linewidth',2);
    %plot(freqbin,smooth(Vscb_Psd,'moving',13)*Vn_Psd(1)/Vscb_Psd(1),'-','Color',lightcoral,'linewidth',2);
    axis([0 100,0,0.1]);
    set(gca,'fontsize',24);
    xlabel('Freq.  [Hz]');
    ylabel('V_{PSD}');
    legend('Current-based','Conductance-based','White noise');

end


PWL = 0;
if PWL == 1

    % Based on 
    %           Response_PWLinearWhiteNoiseInputs.m
    
    % Functions
    
    pwlv=@(v,gL,gc,alpha) v.*(v<alpha)+(alpha+gc/gL*(v-alpha)).*(v>=alpha);
    
    % Parameters
    
    alpha = 1;
    gc = 0.1;
    
    I = 1;
  
    % Eigenvalues of the autonomous system: "broken piece" regime 

    a = -gc/C;
    b = -g1/C;
    c = 1/tau1;
    d = -1/tau1;

    [rc,muc,fnatc] = Eigenvalues2D(a,b,c,d);

    % Impedance profile: "broken piece" regime

    freq = 1:1:100;
    [Zanlc,Phic,fresc,Zmaxc,QZc] = Impedance2D(a,b,c,d,freq);
    
    % Solution to the autonomous system: different initial conditions

    % Time definitions

    Tmax = 1000;
    dt = 0.1;
    t = 0:dt:Tmax;

    va = zeros(1,length(t));
    wa = zeros(1,length(t));

    va(1) = 0;
    wa(1) = 0;
    
    for j=1:length(t)-1
    k1v = -gL*pwlv(va(j),gL,gc,alpha)-g1*wa(j)+I;
    k1w = (va(j)-wa(j))/tau1;
    av = va(j)+k1v*dt;
    aw = wa(j)+k1w*dt;
    k2v = -gL*pwlv(av,gL,gc,alpha)-g1*aw+I;
    k2w = (av-aw)/tau1;
    va(j+1) = va(j)+(k1v+k2v)*dt/2;
    wa(j+1) = wa(j)+(k1w+k2w)*dt/2;
end

    vapwl = va;
    wapwl = wa;

    vv = -100:0.01:100;
    figure(101)
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
    plot(-100,-100,'-g','linewidth',2);
    plot(vv,-gL*pwlv(vv,gL,gc,alpha)/g1,'-r','linewidth',2);
    plot(vv,(-gL*pwlv(vv,gL,gc,alpha)+I)/g1,'--r','linewidth',2);
    plot(vv,vv,'-g','linewidth',2);
    plot(vapwl,wapwl,'-b','linewidth',2);
    plot([-20 20],[0 0],'--k')
    plot([0 0],[-20 20],'--k')
    axis([-2 8 -2 6])
    set(gca,'fontsize',24);
    xlabel('V');
    ylabel('w');
    legend('v-nullclne','w-nullcline','trajectory','Location','NorthEast');  
    
    
    % Poisson inputs
    
    % Time definitions

    Tmax = Tmaxx;
    dt = 0.1;
    t = 0:dt:Tmax;
    
    % Excitatory synaptic-like current-based inputs

    Gsyn = 1;


    % Computation of the numerical solution

    vs = zeros(1,length(t));
    ws = zeros(1,length(t));
     
    for j=1:length(t)-1
        k1v = -gL*pwlv(vs(j),gL,gc,alpha)-g1*ws(j)+Gsyn*Isyntrain(j);
        k1w = (vs(j)-ws(j))/tau1;
        av = vs(j)+k1v*dt;
        av = av+sqrt(2*D*dt)*eta(j);
        aw = ws(j)+k1w*dt;
        k2v = -gL*pwlv(av,gL,gc,alpha)-g1*aw+Gsyn*Isyntrain(j+1);
        k2w = (av-aw)/tau1;
        vs(j+1) = vs(j)+(k1v+k2v)*dt/2;
        vs(j+1) = vs(j+1)+sqrt(2*D*dt)*eta(j);
        ws(j+1) = ws(j)+(k1w+k2w)*dt/2;
    end
    
    vspwl = vs;
    wspwl = ws;

    % Excitatory synaptic-like conductance-based inputs

    Gsyncb = 1;

    % Excitatory 

    Esyn = 1;
    
    % Computation of the numerical solution

    vscb = zeros(1,length(t));
    wscb = zeros(1,length(t));

    for j=1:length(t)-1
        k1v = -gL*pwlv(vscb(j),gL,gc,alpha)-g1*wscb(j)-Gsyncb*Isyntrain(j)*(vscb(j)-Esyn);
        k1w = (vscb(j)-wscb(j))/tau1;
        av = vscb(j)+k1v*dt;
        av = av+sqrt(2*D*dt)*eta(j);
        aw = wscb(j)+k1w*dt;
        k2v = -gL*pwlv(av,gL,gc,alpha)-g1*aw-Gsyncb*Isyntrain(j+1)*(av-Esyn);
        k2w = (av-aw)/tau1;
        vscb(j+1) = vscb(j)+(k1v+k2v)*dt/2;
        vscb(j+1) = vscb(j+1)+sqrt(2*D*dt)*eta(j);
        wscb(j+1) = wscb(j)+(k1w+k2w)*dt/2;
    end

    vscbe = vscb;
    wscbe = wscb;
    
    vscbepwl = vscbe;
    wscbepwl = wscbe;
    
    % Computation of the numerical solution

    v = zeros(1,length(t));
    w = zeros(1,length(t));
    v(1) = 0;
    w(1) = 0;

    for j=1:length(t)-1
        k1v = -gL*pwlv(v(j),gL,gc,alpha)-g1*w(j);
        k1w = (v(j)-w(j))/tau1;
        av = v(j)+k1v*dt;
        av = av+sqrt(2*D*dt)*eta(j);
        aw = w(j)+k1w*dt;
        k2v = -gL*pwlv(av,gL,gc,alpha)-g1*aw;
        k2w = (av-aw)/tau1;
        v(j+1) = v(j)+(k1v+k2v)*dt/2;
        v(j+1) = v(j+1)+sqrt(2*D*dt)*eta(j);
        w(j+1) = w(j)+(k1w+k2w)*dt/2;
    end

    vn = v;
    wn = w;
   
    vnpwl = vn;
    wnpwl = wn;

    % Power spectrum density computation

    [PSD,Freq,PsdManSmooth,freqbin] = powerspectrum(t,Isyntrain);
    ZAPpwl_Psd = PsdManSmooth;

    [PSD,Freq, PsdManSmooth,freqbin] = powerspectrum(t,vspwl);
    Vspwl_Psd = PsdManSmooth;

    Zspwl = Vspwl_Psd./ZAPpwl_Psd;


    [PSD,Freq, PsdManSmooth,freqbin] = powerspectrum(t,vscbepwl);
    Vscbpwl_Psd = PsdManSmooth;

    Zscbpwl = Vscbpwl_Psd./ZAPpwl_Psd;

    [PSD,Freq,PsdManSmooth,freqbin] = powerspectrum(t,eta);
    ZAPnpwl_Psd = PsdManSmooth;

    [PSD,Freq, PsdManSmooth,freqbin] = powerspectrum(t,vnpwl);
    Vnpwl_Psd = PsdManSmooth;


    Znpwl = Vnpwl_Psd./ZAPnpwl_Psd;

    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(-100,-100,'og','linewidth',2);
    plot(freqbin,Zspwl,'ob','linewidth',1);
    plot(freq,Zanl,'-r','linewidth',2);
    plot(freqbin,Znpwl,'-g','linewidth',2);
    axis([0 100,-0.01,5]);
    set(gca,'fontsize',20);
    xlabel('Freq.  [Hz]');
    ylabel('Z');
    legend('Z_{Poisson}','Z_{ANL}','Z_{wn}');



    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(-100,-100,'og','linewidth',2);
    plot(freqbin,Vspwl_Psd,'ob','linewidth',1);
    plot(freqbin,smooth(Vspwl_Psd,'moving',13),'-b','linewidth',2);
    plot(freqbin,Vscbpwl_Psd,'or','linewidth',1);
    plot(freqbin,smooth(Vscbpwl_Psd,'moving',13),'-r','linewidth',2);
    plot(freqbin,Vnpwl_Psd,'og','linewidth',1);
    plot(freqbin,smooth(Vnpwl_Psd,'moving',13),'-g','linewidth',2);
    axis([0 100,-0.01,0.1]);
    set(gca,'fontsize',20);
    xlabel('Freq.  [Hz]');
    ylabel('V_{PSD}');
    legend('Current-based','Conductance-based');



end
