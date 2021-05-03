% 2021-03-20

% Autocorrelations of the amplitudes of the constant pieces in the 
% piecewise constant input functions used in
% Response_LinearWhiteNoiseInputs.m

clearvars;
close all;

Tmax = 100000;
dt = 0.1;
t = 0:dt:Tmax;

Tdur = 1;
Npieces = floor(Tmax/Tdur);
ton = 0:Tdur:Tmax;
ton(1) = dt;
jon = floor(ton/dt);

% Random input amplitudes
    
eta = randn(1,Npieces);
xeta = eta; 
 
% Uniform input amplitudes 

etamax = 2;
etamin = -2;
etaord = etamin:4/Npieces:(Npieces-1)*etamax/Npieces;
P = randperm(length(etaord),length(etaord));
etaperm = etaord(P);
xeta = etaperm;


% Deterministic Gaussian amplitudes

etamax = 2;
etamin = -2;
etaaux = etamin:(etamax-etamin)/(Npieces):(Npieces)*etamax/Npieces;
pd = makedist('Normal','mu',0,'sigma',1);
eta_cdf = cdf(pd,etaaux); 
AmpInt = 2*flip(eta_cdf(1:Npieces/2));
AmpInt = AmpInt*(etamax-etamin)/(2*sum(AmpInt));
etadg(1) = etamin;
for j=2:Npieces/2
    etadg(j) = etadg(j-1)+AmpInt(j-1);
end
for j=1:Npieces/2
    etadg(Npieces/2+j) = -etadg(Npieces/2-j+1);
end
P = randperm(length(etadg),length(etadg));
etaperm = etadg(P);
Xeta= etaperm;
 
Eta_ac = xcorr(xeta,'biased');
 
figure
hold on
plot(Eta_ac,'-b','linewidth',2);
axis([0 2*Npieces+1 -0.1 2]);
set(gca,'fontsize',24);
xlabel('lag');
title('\eta - autocorrelogram');

Ieta = zeros(1,length(t));    
for l=1:length(ton)-1        
    Ieta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = xeta(l);
end
Ieta(end) = Ieta(end-1);

Ieta_ac = xcorr(Ieta,'biased');

figure
hold on
plot(Ieta_ac,'-b','linewidth',2);
axis([0 2*length(Ieta)+1 -0.1 2]);
set(gca,'fontsize',24);
xlabel('lag');
title('I_{\eta} - autocorrelogram');

 
 