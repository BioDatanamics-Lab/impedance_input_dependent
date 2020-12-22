% 2019-08-30

clearvars;
% close all;

%%%input frequency
FF = 1; %[Hz]

vphalf = -38;
vpslp = 6.5;
vrhalf = -79.2;
vrslp =9.78;
C=1;
EL =-65;
ENa =55;
Eh =-20;
gL = 0.5;
gp = 0.1;
gh = 1.5;
Iapp = -2.5; 
tau=80;

IL = @(v) gL*(v - EL);
Ih = @(v) gh*r*(v - Eh);
INap = @(v) gp* pinf(v)*(v - ENa);
pinf = @(v) 1 ./ (1 + exp(-(v-vphalf)./vpslp));
rinf = @(v) 1 ./ (1 + exp((v-vrhalf)./vrslp));

Asinus=2; Asyn=2; Asyn_cond=2; Asqw=2;

Ain = Asinus;
dt = 0.01;
Ncycles = 30;

Eampa=-56.655;

% Generation of the ZAP-like regular and random funtions
Fsin=@(t,f,Ain) Ain*sin(2*pi*f*t/1000-pi/2);
Fsyn=@(t,Ain,t0,taudec) Ain*exp(-(t-t0)/taudec);
Fsqw=@(t,f,Ain,DC) Ain*square(2*pi*f*t/1000-pi/2,DC);

RNDM = 0;
if RNDM == 0
    CycleFreq = 1:1:Ncycles;
elseif RNDM == 1
    CycleFreq = randperm(Ncycles,Ncycles);
end

CycleFreq = FF*ones(1,length(CycleFreq));

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
v(1)=-58;
w(1)=0.1;

for j=1:length(t)-1
    k1v = -gL*(v(j)-EL) -gh*w(j)*(v(j)-Eh) -gp.*pinf(v(j)).*(v(j)-ENa) + Izaplike(j) + Iapp; 
    k1w = (rinf(v(j)) - w(j))/tau; 
    av = v(j)+k1v*dt;
    aw = w(j)+k1w*dt;
    k2v = -gL*(av-EL) -gh*aw*(av-Eh) -gp.*pinf(av).*(av-ENa) + Izaplike(j+1) + Iapp; 
    k2w = (rinf(av) - aw)/tau; 
    v(j+1) = v(j)+(k1v+k2v)*dt/2;
    w(j+1) = w(j)+(k1w+k2w)*dt/2;
end

% Generation of the ZAP-like regular and random synaptic-like funtions
Ain = Asyn;
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
vs(1)=-58;
ws(1)=0.1;

for j=1:length(t)-1
    k1v = -gL*(vs(j)-EL) -gh*ws(j)*(vs(j)-Eh) -gp.*pinf(vs(j)).*(vs(j)-ENa) + Isyntrain(j) + Iapp; 
    k1w = (rinf(vs(j)) - ws(j))/tau; 
    av = vs(j)+k1v*dt;
    aw = ws(j)+k1w*dt;
    k2v = -gL*(av-EL) -gh*aw*(av-Eh) -gp.*pinf(av).*(av-ENa) + Isyntrain(j+1) + Iapp; 
    k2w = (rinf(av) - aw)/tau; 
    vs(j+1) = vs(j)+(k1v+k2v)*dt/2;
    ws(j+1) = ws(j)+(k1w+k2w)*dt/2;
end


%%% conductance-based

Ain = Asyn_cond;
taurse_e = 1;     % AMPA
taudec_e = 5;     % AMPA
taurse_i = 0.2;   % GABA_A     
taudec_i = 10;    % GABA_A
Eex = 20;
Ein = -60;

taudec = taudec_e;

t = 0:dt:CyclePer(1);
Isyntrainc = Fsyn(t,Ain,t(1),taudec);
for k=2:Ncycles
    tloc = dt:dt:CyclePer(k);  
    Floc = Fsyn(tloc,Ain,tloc(1),taudec);
    t = [t t(end)+tloc];
    Isyntrainc = [Isyntrainc Floc];
end

% Computation of the numerical solution

vsc = zeros(1,length(t));
wsc = zeros(1,length(t));
vsc(1)=-58;
wsc(1)=0.1;
condinput = zeros(1,length(t));

for j=1:length(t)-1
    k1v = -gL*(vsc(j)-EL) -gh*wsc(j)*(vsc(j)-Eh) -gp.*pinf(vsc(j)).*(vsc(j)-ENa) + Isyntrainc(j)*(Eampa-vsc(j)) + Iapp; 
    k1w = (rinf(vsc(j)) - wsc(j))/tau; 
    av = vsc(j)+k1v*dt;
    aw = wsc(j)+k1w*dt;
    k2v = -gL*(av-EL) -gh*aw*(av-Eh) -gp.*pinf(av).*(av-ENa) + Isyntrainc(j+1)*(Eampa-av) + Iapp; 
    k2w = (rinf(av) - aw)/tau; 
    vsc(j+1) = vsc(j)+(k1v+k2v)*dt/2;
    wsc(j+1) = wsc(j)+(k1w+k2w)*dt/2;
    condinput(j) = Isyntrainc(j)*(Eampa-vsc(j));
end


% Generation of the ZAP-like square-wave funtions
Ain = Asqw;
DC = 50;

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
vsqw(1)=-58;
wsqw(1)=0.1;

for j=1:length(t)-1
    k1v = -gL*(vsqw(j)-EL) -gh*wsqw(j)*(vsqw(j)-Eh) -gp.*pinf(vsqw(j)).*(vsqw(j)-ENa) + Isqwtrain(j) + Iapp; 
    k1w = (rinf(vsqw(j)) - wsqw(j))/tau; 
    av = vsqw(j)+k1v*dt;
    aw = wsqw(j)+k1w*dt;
    k2v = -gL*(av-EL) -gh*aw*(av-Eh) -gp.*pinf(av).*(av-ENa) + Isqwtrain(j+1) + Iapp; 
    k2w = (rinf(av) - aw)/tau; 
    vsqw(j+1) = vsqw(j)+(k1v+k2v)*dt/2;
    wsqw(j+1) = wsqw(j)+(k1w+k2w)*dt/2;
end


figure; 
subplot(4,2,1)
plot(t,Izaplike,'-b','linewidth',2);
xlim([0 Tmax]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('input');

subplot(4,2,2)
plot(t,v,'-b','linewidth',2);
xlim([0 Tmax]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('v [mV]');

subplot(4,2,3)
plot(t,Isyntrain,'-b','linewidth',2);
xlim([0 Tmax])
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('input');

subplot(4,2,4)
plot(t,vs,'-b','linewidth',2);
xlim([0 Tmax])
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('v [mV]');

subplot(4,2,5)
plot(t,condinput,'-b','linewidth',2);
xlim([0 Tmax])
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('input');

subplot(4,2,6)
plot(t,vsc,'-b','linewidth',2);
xlim([0 Tmax])
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('v [mV]');

subplot(4,2,7)
plot(t,Isqwtrain,'-b','linewidth',2);
xlim([0 Tmax])
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('input');

subplot(4,2,8)
plot(t,vsqw,'-b','linewidth',2);
xlim([0 Tmax])
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('v [mV]');