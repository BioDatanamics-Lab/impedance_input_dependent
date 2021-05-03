% 2021-03-26

% Response of PWL systems to white noise and PWC noise
%
% Based on Response_LinearWhiteNoiseInputs.m (PWC = 3)

% Based

clearvars;
close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightsalmon = [0.9 0.5 0.4];
mediumacquamarine = [0.4 0.8 0.6];
lightgray = [.7 .7 .7];
darkgray = [.3 .3 .3];
darkgray2 = [.1 .1 .1];

% Functions

heaviside=@(t) 0.5*(t == 0)+(t > 0);
pwlv=@(v,gL,gc,alpha) v.*(v<alpha)+(alpha+gc/gL*(v-alpha)).*(v>=alpha);
% gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
gaus = @(x,mu,sig) 1/(sig*sqrt(2*pi))*exp(-(((x-mu).^2)/(2*sig.^2)));


C = 1;
gL = 0.25;
g1 = 0.25;
tau1 = 100;

alpha = 1;
gc = 0.1;

Gsyn = 2;
Esyn = 1;

I = 1;

D = 0.4;

% Eigenvalues of the autonomous system: piece 1

a = -gL/C;
b = -g1/C;
c = 1/tau1;
d = -1/tau1;

[r,mu,fnat] = Eigenvalues2D(a,b,c,d);

% Impedance profile: piece 1

freq = 1:1:100;
[Zanl,Phi,fres,Zmax,QZ] = Impedance2D(a,b,c,d,freq);

% Eigenvalues of the autonomous system: piece 1

a = -gc/C;
b = -g1/C;
c = 1/tau1;
d = -1/tau1;

[rc,muc,fnatc] = Eigenvalues2D(a,b,c,d);

% Impedance profile: piece 1

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
    k1v = -gL*va(j)-g1*wa(j)+I;
    k1w = (va(j)-wa(j))/tau1;
    av = va(j)+k1v*dt;
    aw = wa(j)+k1w*dt;
    k2v = -gL*av-g1*aw+I;
    k2w = (av-aw)/tau1;
    va(j+1) = va(j)+(k1v+k2v)*dt/2;
    wa(j+1) = wa(j)+(k1w+k2w)*dt/2;
end

valin = va;
walin = wa;

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

for j=1:length(t)-1
    k1v = -gL*pwlv(va(j),gL,gc,alpha)-g1*wa(j)-Gsyn*I*(va(j)-Esyn);
    k1w = (va(j)-wa(j))/tau1;
    av = va(j)+k1v*dt;
    aw = wa(j)+k1w*dt;
    k2v = -gL*pwlv(av,gL,gc,alpha)-g1*aw-Gsyn*I*(av-Esyn);
    k2w = (av-aw)/tau1;
    va(j+1) = va(j)+(k1v+k2v)*dt/2;
    wa(j+1) = wa(j)+(k1w+k2w)*dt/2;
end

vacbpwl = va;
wacbpwl = wa;

figure
hold on
plot(t,vapwl,'-b','linewidth',2);
plot(t,valin,'--b','linewidth',2);
plot(t,vacbpwl,'-','Color',lightblueish,'linewidth',2);
axis([0 200 0 2]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('V');


vv = -100:0.01:100;
figure
hold on
plot(vv,-gL*vv/g1,'-r','linewidth',2);
plot(vv,vv,'-g','linewidth',2);
plot(valin,walin,'-b','linewidth',2);
plot(vv,(-gL*vv+I)/g1,'--r','linewidth',2);
plot([-20 20],[0 0],'--k')
plot([0 0],[-20 20],'--k')
axis([-2 8 -2 6])
set(gca,'fontsize',24);
xlabel('V  [mV]');
ylabel('w');
%legend('v-nullclne','w-nullcline','trajectory','Location','NorthEast');

figure
hold on
plot(vv,-gL*pwlv(vv,gL,gc,alpha)/g1,'-r','linewidth',2);
plot(vv,vv,'-g','linewidth',2);
plot(vapwl,wapwl,'-b','linewidth',2);
plot(vv,(-gL*pwlv(vv,gL,gc,alpha)+I)/g1,'--r','linewidth',2);
plot([-20 20],[0 0],'--k')
plot([0 0],[-20 20],'--k')
axis([-2 8 -2 6])
set(gca,'fontsize',24);
xlabel('V  [mV]');
ylabel('w');
%legend('v-nullclne','w-nullcline','trajectory','Location','NorthEast');

figure
hold on
plot(vv,-gL*pwlv(vv,gL,gc,alpha)/g1,'-r','linewidth',2);
plot(vv,vv,'-g','linewidth',2);
plot(vacbpwl,wacbpwl,'-b','linewidth',2);
plot(vv,(-gL*pwlv(vv,gL,gc,alpha)-Gsyn*I*(vv-Esyn))/g1,'--r','linewidth',2);
plot([-20 20],[0 0],'--k')
plot([0 0],[-20 20],'--k')
axis([-2 8 -2 6]);
set(gca,'fontsize',24);
xlabel('V  [mV]');
ylabel('w');
%legend('v-nullclne','w-nullcline','trajectory','Location','NorthEast');


Tdur = 1;
Npieces = floor(Tmax/Tdur);
ton = 0:Tdur:Tmax;
ton(1) = dt;
jon = floor(ton/dt);


% Uniform input amplitudes

etamax = 2;
etamin = 0;
eta = etamin:(etamax-etamin)/Npieces:(Npieces-1)*etamax/Npieces;
eta = eta(randperm(length(eta),length(eta)));

Ieta = zeros(1,length(t));    
for l=1:length(ton)-1        
    Ieta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = eta(l);
end
Ieta(end) = Ieta(end-1);

v = zeros(1,length(t));
w = zeros(1,length(t));

v(1) = mean(eta);
w(1) = mean(eta);


for j=1:length(t)-1
    k1v = -gL*v(j)-g1*w(j)+D*Ieta(j);
    k1w = (v(j)-w(j))/tau1;
    av = v(j)+k1v*dt;
    aw = w(j)+k1w*dt;
    k2v = -gL*av-g1*aw+D*Ieta(j+1);
    k2w = (av-aw)/tau1;
    v(j+1) = v(j)+(k1v+k2v)*dt/2;
    w(j+1) = w(j)+(k1w+k2w)*dt/2;
end

vetalin = v;
wetalin = w;

figure
hold on
plot(t,vetalin,'-b','linewidth',2);
plot(t,0.25*Ieta-4,'-','Color',lightgray,'linewidth',1);
axis([0 Tmax -5 5]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('V  [mV]');

for j=1:length(t)-1
    k1v = -gL*pwlv(v(j),gL,gc,alpha)-g1*w(j)+D*Ieta(j);
    k1w = (v(j)-w(j))/tau1;
    av = v(j)+k1v*dt;
    aw = w(j)+k1w*dt;
    k2v = -gL*pwlv(av,gL,gc,alpha)-g1*aw+D*Ieta(j+1);
    k2w = (av-aw)/tau1;
    v(j+1) = v(j)+(k1v+k2v)*dt/2;
    w(j+1) = w(j)+(k1w+k2w)*dt/2;
end

vetapwl = v;
wetapwl = w;

figure
hold on
plot(t,vetapwl,'-b','linewidth',2);
plot(t,0.25*Ieta-4,'-','Color',lightgray,'linewidth',1);
axis([0 Tmax -5 5]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('V  [mV]');

for j=1:length(t)-1
    k1v = -gL*pwlv(v(j),gL,gc,alpha)-g1*w(j)-Gsyn*D*Ieta(j)*(v(j)-Esyn);
    k1w = (v(j)-w(j))/tau1;
    av = v(j)+k1v*dt;
    aw = w(j)+k1w*dt;
    k2v = -gL*pwlv(av,gL,gc,alpha)-g1*aw-Gsyn*D*Ieta(j+1)*(av-Esyn);
    k2w = (av-aw)/tau1;
    v(j+1) = v(j)+(k1v+k2v)*dt/2;
    w(j+1) = w(j)+(k1w+k2w)*dt/2;
end

vetacbpwl = v;
wetacbpwl = w;

figure
hold on
plot(t,vetacbpwl,'-b','linewidth',2);
plot(t,0.25*Ieta-4,'-','Color',lightgray,'linewidth',1);
axis([0 Tmax -5 5]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('V  [mV]');


PWS = 0;
if PWS == 1
    
   % Time definitions

    Tmaxx = 1000000;
    dt = 0.1;
    t = 0:dt:Tmaxx;

    Tdur = 1;
    Npieces = floor(Tmaxx/Tdur);
    ton = 0:Tdur:Tmaxx;
    ton(1) = dt;
    jon = floor(ton/dt);


    % Uniform input amplitudes

    etamax = 2;
    etamin = 0;
    eta = etamin:(etamax-etamin)/Npieces:(Npieces-1)*etamax/Npieces;
    eta = eta(randperm(length(eta),length(eta)));
    

    Ieta = zeros(1,length(t));    
    for l=1:length(ton)-1        
        Ieta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = eta(l);
    end
    Ieta(end) = Ieta(end-1);
    
    v(1) = mean(eta);
    w(1) = mean(eta);


    for j=1:length(t)-1
        k1v = -gL*v(j)-g1*w(j)+D*Ieta(j);
        k1w = (v(j)-w(j))/tau1;
        av = v(j)+k1v*dt;
        aw = w(j)+k1w*dt;
        k2v = -gL*av-g1*aw+D*Ieta(j+1);
        k2w = (av-aw)/tau1;
        v(j+1) = v(j)+(k1v+k2v)*dt/2;
        w(j+1) = w(j)+(k1w+k2w)*dt/2;
    end

    vetalin = v;
    wetalin = w;
    
    [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,D*Ieta);
    ZAP_Psd = PsdManSmooth;

    [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vetalin);
    V_Psd = PsdManSmooth;

    Z = V_Psd./ZAP_Psd;
    
    Vlin_Psd = V_Psd;
    Zlin = Z;

    for j=1:length(t)-1
        k1v = -gL*pwlv(v(j),gL,gc,alpha)-g1*w(j)+D*Ieta(j);
        k1w = (v(j)-w(j))/tau1;
        av = v(j)+k1v*dt;
        aw = w(j)+k1w*dt;
        k2v = -gL*pwlv(av,gL,gc,alpha)-g1*aw+D*Ieta(j+1);
        k2w = (av-aw)/tau1;
        v(j+1) = v(j)+(k1v+k2v)*dt/2;
        w(j+1) = w(j)+(k1w+k2w)*dt/2;
    end

    vetapwl = v;
    wetapwl = w;
    
    [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,D*Ieta);
    ZAP_Psd = PsdManSmooth;

    [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vetapwl);
    V_Psd = PsdManSmooth;

    Z = V_Psd./ZAP_Psd;
    
    Vpwl_Psd = V_Psd;
    Zpwl = Z;

 
    for j=1:length(t)-1
        k1v = -gL*pwlv(v(j),gL,gc,alpha)-g1*w(j)-Gsyn*D*Ieta(j)*(v(j)-Esyn);
        k1w = (v(j)-w(j))/tau1;
        av = v(j)+k1v*dt;
        aw = w(j)+k1w*dt;
        k2v = -gL*pwlv(av,gL,gc,alpha)-g1*aw-Gsyn*D*Ieta(j+1)*(av-Esyn);
        k2w = (av-aw)/tau1;
        v(j+1) = v(j)+(k1v+k2v)*dt/2;
        w(j+1) = w(j)+(k1w+k2w)*dt/2;
    end

    vetacbpwl = v;
    wetacbpwl = w;
    
    [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,D*Ieta);
    ZAP_Psd = PsdManSmooth;

    [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vetacbpwl);
    V_Psd = PsdManSmooth;

    Z = V_Psd./ZAP_Psd;
    
    Vcbpwl_Psd = V_Psd;
    Zcbpwl = Z;
    
    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(freqbin,Zlin,'ob','linewidth',1);
    plot(freqbin,Vlin_Psd*250,'or','linewidth',1);
%     plot(fres,0,'ok','linewidth',3);
%     plot(fnat,0,'ok','linewidth',3);
%     plot([fres fres],[0 20],'--k','linewidth',2);
%     plot([fnat fnat],[0 20],'--k','linewidth',2);
    %plot(freqbin,Z,'-b','linewidth',1);
    axis([0 100,-0.01,6]);
    set(gca,'fontsize',20);
    xlabel('Freq.  [Hz]');
    ylabel('Z');
    
    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(freqbin,Zpwl,'ob','linewidth',1);
    plot(freqbin,Vpwl_Psd*250,'or','linewidth',1);
%     plot(fres,0,'ok','linewidth',3);
%     plot(fnat,0,'ok','linewidth',3);
%     plot([fres fres],[0 20],'--k','linewidth',2);
%     plot([fnat fnat],[0 20],'--k','linewidth',2);
    %plot(freqbin,Z,'-b','linewidth',1);
    axis([0 100,-0.01,6]);
    set(gca,'fontsize',20);
    xlabel('Freq.  [Hz]');
    ylabel('Z');
    
    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(freqbin,Zcbpwl,'ob','linewidth',1);
    plot(freqbin,Vcbpwl_Psd*250,'or','linewidth',1);
%     plot(fres,0,'ok','linewidth',3);
%     plot(fnat,0,'ok','linewidth',3);
%     plot([fres fres],[0 20],'--k','linewidth',2);
%     plot([fnat fnat],[0 20],'--k','linewidth',2);
    %plot(freqbin,Z,'-b','linewidth',1);
    axis([0 100,-0.01,6]);
    set(gca,'fontsize',20);
    xlabel('Freq.  [Hz]');
    ylabel('Z');
    
    figure
    hold on
    plot(freqbin,Vpwl_Psd,'ob','linewidth',1);
    plot(freqbin,Vlin_Psd,'or','linewidth',1);
    plot(freqbin,Vcbpwl_Psd,'og','linewidth',1);
    axis([0 100,-0.01,6]);
    set(gca,'fontsize',20);
    xlabel('Freq.  [Hz]');
    ylabel('V_Psd');
    
    
end




