% 2021-01-27
% 2020-10-16
% 2019-08-30
 
% Cell responses to Zap-like sinusoidal (SIN), square-wave (SQW) and
% synaptic-like (SYN) inputs (excitatory, inhibitory, current- and
% conductance-based). 
% Each signal consist of a number of cycles with different frequencies
% regularly ordered (increasing order of frequencies with the cycle number) 
% or randomly ordered. 
% The response frequency can be comnputed for each input cycle (using the
% input cycles as reference) or for each output cycle (peak responses
% belong to frequencies computed as normalized differences between the
% adjacent trough times and trough responses belong to frequencies computed
% as normalized differences between the adjacent peak times.
% 
% 

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

Fsin=@(t,f,Ain) Ain*sin(2*pi*f*t/1000-pi/2);
Fsyn=@(t,Ain,t0,taudec) Ain*exp(-(t-t0)/taudec);
Fsqw=@(t,f,Ain,DC) Ain*square(2*pi*f*t/1000-pi/2,DC);


C = 1;
gL = 0.25;
g1 = 0.25;
tau1 = 100;

% gL = 0.05;
% g1 = 0.3;

% g1 = 0;
%

% C = 1;
% gL = 0.1;
% g1 = 0.8;
% tau1 = 100;


D = 0;

a = -gL/C;
b = -g1/C;
c = 1/tau1;
d = -1/tau1;


Ain = 1;
dt = 0.01;
Ncycles = 100;


% Eigenvalues of the autonomous system 

[r,mu,fnat] = Eigenvalues2D(a,b,c,d);

% Impedance profile

[Z,Phi,fres,Zmax,QZ] = Impedance2D(a,b,c,d,[1:1:100]);


% Generation of the ZAP-like regular and random cycle periods

RNDM = 0;
if RNDM == 0
    CycleFreq = 1:1:Ncycles;
elseif RNDM == 1
    CycleFreq = randperm(Ncycles,Ncycles);
end
CyclePer = 1000./CycleFreq;



% Type of input 
%               SIN: sinusoidal
%               SQW: square-wave, 
%               SYN: synaptic-like  
%                    1: single trial (LIN)
%                    2: multiple trials (LIN)
%                    3: single trial (PWL)
            
SIN = 0;
SQW = 0;
SYN = 3; 


% Ouput cycle period (local frequency) computation


            % NOTE HGR. There are two
            % types of oscillatory peaks/troughs: (1) the envelope ones,
            % and (2) the low-amplitude ones that arise in the response
            % mixed-mode oscillatory (MMO) patterns and result from the
            % transient dynamics in response to the oscillatory inputs.
            % They are seen in response to square-wave and synaptic-like
            % patterns (and less or none in resposne to sinusoidal inputs).
            % There are various ways of computing the frequency-ampliitude
            % response profiles, which depend on the different cases.
            % (1) Sometimes it is possible to use the input cycles as a
            % reference and (2) others they need to be referred to the
            % input cycles. The differences between the two are mostly 
            % technical when there are no MMO responses. When there are MMO
            % responses, one needs to use the two simultaneously. The large
            % amplitude peaks/troughs will be captured twice. This will not
            % affect the graphs, but it might affect the quantification if
            % the peak/trough data is used for that purpose. To be more
            % specific: 
            % (1) Within the boundaries of each input cycle, compute the
            % max and min of the response.
            % (2) Along (all the points of the response, compute the max 
            % and min.
            % Once the peaks and troughs have been identified, the need to
            % be ordered according to increasing values of the frequency.
            % To this end we use a number of functions:
            % (i) OrderingCycleFreq
            % (ii) OrderingCycleFreqPeak
            % (iii) OrderingCycleFreqTrough
            % (iv) OrderingCycleFreqPeakTrough
            % If possible and meaningful, we also compute the phase. 
            %
          

            
if SIN == 1

    %%%%% Generation of the ZAP-like sinusoidl funtions

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
    Deltajbdry = diff(jbdry);
    
    % White (Gaussian) noise

    eta = randn(1,length(t));

    % Computation of the numerical solution

    v = zeros(1,length(t));
    w = zeros(1,length(t));
    v(1) = -Ain/(gL+g1);
    w(1) = v(1);


    for j=1:length(t)-1
        k1v = -gL*v(j)-g1*w(j)+Izaplike(j);
        k1w = (v(j)-w(j))/tau1;
        av = v(j)+k1v*dt;
        av = av+sqrt(2*D*dt)*eta(j);
        aw = w(j)+k1w*dt;
        k2v = -gL*av-g1*aw+Izaplike(j+1);
        k2w = (av-aw)/tau1;
        v(j+1) = v(j)+(k1v+k2v)*dt/2;
        v(j+1) = v(j+1)+sqrt(2*D*dt)*eta(j);
        w(j+1) = w(j)+(k1w+k2w)*dt/2;
    end

    % Computation of the peaks and peak-times
    % no trough and trough times computation (vmin-->vmax)

    
    [tpeakzap,zapmax,contpeakzap] = PeaksOsc(Izaplike,t,dt,Tmax,dt);
    [tpeak,vmax,cntpeak] = PeaksOsc(v,t,dt,Tmax,dt);
    [ttrough,vmin,cnttrough] = TroughsOsc(v,t,dt,Tmax,dt);
    ttrough = [0 ttrough];
    vmin = [v(1) vmin];
    Phase = zeros(1,Ncycles);
    for k=1:Ncycles
        Phase(k) = (tpeak(k)-tpeakzap(k))*CycleFreq(k);
    end
    Phase = Phase*2*pi/1000;
    CycleFreqOutputPeak = 1000./diff(ttrough);
    [CycleFreqOrdPeak,vmaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vmax(1:end-1));
    CycleFreqOutputTrough = 1000./diff(tpeak);
    [CycleFreqOrdTrough,vminord] = OrderingCycleFreqTrough(CycleFreqOutputTrough,vmin(2:end));
        
   

    vsinmaxord = vmaxord;
    vsinminord = vminord;
    tsinpeak = tpeak;
    tsintrough = ttrough; 

    figure
    hold on
    plot(t,v,'-b','linewidth',2);
    plot(t,Izaplike-8,'-','Color',[.7 .7 .7],'linewidth',2);
    plot(tpeak,vmax,'or')
    plot(ttrough,vmin,'or');
    axis([0 Tmax -10 5]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('');

    
    
%         figure
%         hold on
%         plot(CycleFreqOrd,vmaxord,'ob','linewidth',2);
%         plot(CycleFreqOrd,vminord,'ob','linewidth',2);
%         %plot(CycleFreqOrd,(vmaxord-vminord)/2,'or','linewidth',2);
%         axis([0 100 -5 5]);
%         set(gca,'fontsize',24);
%         xlabel('f  [Hz]');
%         ylabel('v_{max}  [mV]');

    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(CycleFreqOrdPeak,vmaxord,'ob','linewidth',2);
    plot(CycleFreqOrdTrough,vminord,'ob','linewidth',2);
    %plot(CycleFreqOrd,(vmaxord-vminord)/2,'or','linewidth',2);
    axis([0 100 -5 5]);
    set(gca,'fontsize',24);
    xlabel('f  [Hz]');
    ylabel('V_{ENV}^{+/-}  [mV]');
        
    
%     figure
%     hold on
%     plot(CycleFreqOrdPeak,'ob','linewidth',2);
%     plot(CycleFreqOrdTrough,'or','linewidth',2);
%     set(gca,'fontsize',24);
%     xlabel('Cycle #')
%     ylabel('Output f  [Hz]');
    
    % figure
    % hold on
    % plot(t,Izaplike,'-b','linewidth',2);
    % axis([0 Tmax -1.2 1.2]);
    % set(gca,'fontsize',24);
    % xlabel('t');
    % ylabel('');

    % ff = CycleFreqOrd;
    % figure
    % hold on
    % plot(CycleFreqOrd,Phi,'-b','linewidth',2);
    % plot(CycleFreqOrd,Phase,'ob','linewidth',2);
    % plot([0 max(ff)],[0 0],'--','Color',[.7 .7 .7]);
    % plot([0 max(ff)],[pi/2 pi/2],'--','Color',[.7 .7 .7]);
    % plot([0 max(ff)],[-pi/2 -pi/2],'--','Color',[.7 .7 .7]);
    % plot([0 max(ff)],[2*pi 2*pi],'--','Color',[.7 .7 .7]);
    % axis([0 ff(end) -pi/2-0.2 pi/2+0.2])
    % set(gca,'fontsize',24);
    % xlabel('f  [Hz]');
    % ylabel('\Phi');
    % title('Phase profile');
  
end

if SQW == 1

    % Generation of the ZAP-like square-wave funtions

    DC = 50;


    t = 0:dt:CyclePer(1);
    Isqwtrain = Fsqw(t,CycleFreq(1),Ain,DC);
    tbdry(1) = 0;
    jbdry(1) = 1;
    tbdry(2) = t(end);
    jbdry(2) = length(t);
    for k=2:Ncycles
        tloc = dt:dt:CyclePer(k);  
        Floc = Fsqw(tloc,CycleFreq(k),Ain,DC);   
        t = [t t(end)+tloc];
        Isqwtrain = [Isqwtrain Floc];
        tbdry(k+1) = t(end);
        jbdry(k+1) = length(t);
    end
    Tmax = t(end);
    Deltajbdry = diff(jbdry);
    
    % White (Gaussian) noise

    eta = randn(1,length(t));

    % Computation of the numerical solution

    vsqw = zeros(1,length(t));
    wsqw = zeros(1,length(t));
    vsqw(1) = -Ain/(gL+g1);
    wsqw(1) = vsqw(1);


    for j=1:length(t)-1
        k1v = -gL*vsqw(j)-g1*wsqw(j)+Isqwtrain(j);
        k1w = (vsqw(j)-wsqw(j))/tau1;
        av = vsqw(j)+k1v*dt;
        av = av+sqrt(2*D*dt)*eta(j);
        aw = wsqw(j)+k1w*dt;
        k2v = -gL*av-g1*aw+Isqwtrain(j+1);
        k2w = (av-aw)/tau1;
        vsqw(j+1) = vsqw(j)+(k1v+k2v)*dt/2;
        vsqw(j+1) = vsqw(j+1)+sqrt(2*D*dt)*eta(j);
        wsqw(j+1) = wsqw(j)+(k1w+k2w)*dt/2;
    end

    % Computation of the peaks and peak-times

  
    tpeakzap = tbdry(1:Ncycles)+1000./CycleFreq*0.25;
    [tpeak,vsqwmax,cntpeak] = PeaksOsc(vsqw,t,dt,Tmax,dt);
    [ttrough,vsqwmin,cnttrough] = TroughsOsc(vsqw,t,dt,Tmax,dt);
    ttrough = [0 ttrough];
    vsqwmin = [vsqw(1) vsqwmin];
    Phase = zeros(1,Ncycles);
    for k=1:Ncycles
        Phase(k) = (tpeak(k)-tpeakzap(k))*CycleFreq(k);
    end
    Phase = Phase*2*pi/1000;
    CycleFreqOutputPeak = 1000./diff(ttrough);
    [CycleFreqOrdPeak,vmaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsqwmax(1:end-1));
    CycleFreqOutputTrough = 1000./diff(tpeak);
    [CycleFreqOrdTrough,vminord] = OrderingCycleFreqTrough(CycleFreqOutputTrough,vsqwmin(2:end));
    vsqwmax2 = zeros(1,Ncycles);
    vsqwmin2 = zeros(1,Ncycles);
    tsqwpeak2 = zeros(1,Ncycles);
    tsqwtrough2 = zeros(1,Ncycles);
    for k=1:Ncycles
        [vsqwmax2(k),jpeak2(k)] = max(vsqw(jbdry(k):jbdry(k+1)));  
        [vsqwmin2(k),jtrough2(k)] = min(vsqw(jbdry(k):jbdry(k+1))); 
        tsqwpeak2(k) = tbdry(k)+jpeak2(k)*dt;
        tsqwtrough2(k) = tbdry(k)+jtrough2(k)*dt;
    end
    [CycleFreqOrd2,vsqwmaxord2,vsqwminord2] = OrderingCycleFreqPeakTrough(CycleFreq,vsqwmax2,vsqwmin2);
   

    vsqwmaxord = vmaxord;
    vsqwminord = vminord;
    tsqwpeak = tpeak;
    tsqwtrough = ttrough;

    % [peaktzap,zapmax,contpeakzap] = PeaksOsc(Isqwtrain,t,dt,Tmax,dt);
    % [peakt,vsqwmax,cntpeak] = PeaksOsc(vsqw,t,dt,Tmax,dt);
    % [trought,vsqwmin,cnttrough] = TroughsOsc(vsqw,t,dt,Tmax,dt);
    % Phasesqw = zeros(1,Ncycles);
    % for k=1:Ncycles
    %     Phase(k) = (peakt(k)-peaktzap(k))*CycleFreq(k);
    % end
    % Phase = Phase*2*pi/1000;



    % figure
    % hold on
    % plot(t,Isqwtrain,'-b','linewidth',2);
    % axis([0 Tmax -1.2 1.2]);
    % set(gca,'fontsize',24);
    % xlabel('t');
    % ylabel('');

    figure
    hold on
    plot(t,vsqw,'-b','linewidth',2);
    plot(t,Isqwtrain-8,'-','Color',[.7 .7 .7],'linewidth',2);
    plot(tpeak,vsqwmax,'or')
    plot(ttrough,vsqwmin,'or')
    axis([0 Tmax -10 6]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('');

    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(CycleFreqOrdPeak,vsqwmaxord,'ob','linewidth',2);
    plot(CycleFreqOrdTrough,vsqwminord,'ob','linewidth',2);
%     plot(CycleFreqOrd2,vsqwmaxord2,'or','linewidth',2);
%     plot(CycleFreqOrd2,vsqwminord2,'or','linewidth',2);
    %plot(CycleFreqOrd,(vsqwmaxord-vsqwminord)/2,'or','linewidth',2);
    axis([0 100 -10 10]);
    set(gca,'fontsize',24);
    xlabel('f  [Hz]');
    ylabel('v_{max}  [mV]');

end

if SYN == 1

    % Generation of the ZAP-like regular and random synaptic-like funtions

    taurse_e = 1;     % AMPA
    taudec_e = 5;     % AMPA
    taurse_i = 0.2;   % GABA_A     
    taudec_i = 10;    % GABA_A
    Eex = 20;
    Ein = -60;

    taudec = taudec_e;
    

    t = 0:dt:CyclePer(1);
    Isyntrain = Fsyn(t,Ain,t(1),taudec);
    tbdry(1) = 0;
    jbdry(1) = 1;
    tbdry(2) = t(end);
    jbdry(2) = length(t);
    for k=2:Ncycles
        tloc = dt:dt:CyclePer(k);  
        Floc = Fsyn(tloc,Ain,tloc(1),taudec);
        t = [t t(end)+tloc];
        Isyntrain = [Isyntrain Floc];
        tbdry(k+1) = t(end);
        jbdry(k+1) = length(t);
    end
    Tmax = t(end);
    Deltajbdry = diff(jbdry);
    
    % White (Gaussian) noise

    eta = randn(1,length(t));

    % Excitatory synaptic-like current-based inputs

    Gsyn = 1;

    % Computation of the numerical solution

    vs = zeros(1,length(t));
    ws = zeros(1,length(t));

    for j=1:length(t)-1
        k1v = -gL*vs(j)-g1*ws(j)+Gsyn*Isyntrain(j);
        k1w = (vs(j)-ws(j))/tau1;
        av = vs(j)+k1v*dt;
        av = av+sqrt(2*D*dt)*eta(j);
        aw = ws(j)+k1w*dt;
        k2v = -gL*av-g1*aw+Gsyn*Isyntrain(j+1);
        k2w = (av-aw)/tau1;
        vs(j+1) = vs(j)+(k1v+k2v)*dt/2;
        vs(j+1) = vs(j+1)+sqrt(2*D*dt)*eta(j);
        ws(j+1) = ws(j)+(k1w+k2w)*dt/2;
    end


    % Computation of the peaks and peak-times

    [tpeakzap,zapmax,contpeakzap] = PeaksOsc(Isyntrain,t,dt,Tmax,dt);
    tpeakzap = [0 tpeakzap];
    zapmax = [1 zapmax];
    [tsepeak,vsemax,cntpeak] = PeaksOsc(vs,t,dt,Tmax,dt);
    [tsetrough,vsemin,cnttrough] = TroughsOsc(vs,t,dt,Tmax,dt);
    tsetrough = [tsetrough t(end)];
    vsemin = [vsemin vs(end)];
    Phases = zeros(1,Ncycles);
    for k=1:Ncycles
        Phases(k) = (tsepeak(k)-tpeakzap(k))*CycleFreq(k);
    end
    vsemax2 = zeros(1,Ncycles);
    vsemin2 = zeros(1,Ncycles);
    tsepeak2 = zeros(1,Ncycles);
    tsetrough2 = zeros(1,Ncycles);
    for k=1:Ncycles
        [vsemax2(k),jpeak2(k)] = max(vs(jbdry(k):jbdry(k+1)-1));  
        [vsemin2(k),jtrough2(k)] = min(vs(jbdry(k):jbdry(k+1)-1)); 
        tsepeak2(k) = tbdry(k)+jpeak2(k)*dt;
        tsetrough2(k) = tbdry(k)+jtrough2(k)*dt;
    end
    [CycleFreqOrd2,vsemaxord2,vseminord2] = OrderingCycleFreqPeakTrough(CycleFreq,vsemax2,vsemin2);
    CycleFreqOutputPeak = 1000./diff(tsetrough);
    if length(CycleFreqOutputPeak) < length(vsemax)   
        [CycleFreqOrdPeak,vsemaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsemax(1:end-1));
    else
        [CycleFreqOrdPeak,vsemaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsemax);
    end
    CycleFreqOutputTrough = 1000./diff(tsepeak);
    [CycleFreqOrdTrough,vseminord] = OrderingCycleFreqTrough(CycleFreqOutputTrough,vsemin(2:end));     
 

    vse = vs;
    
    

    % figure
    % hold on
    % plot(t,Isyntrain,'-b','linewidth',2);
    % axis([0 Tmax -1.2 1.2]);
    % set(gca,'fontsize',24);
    % xlabel('t');
    % ylabel('');

    figure
    hold on
    plot(t,vse,'-b','linewidth',2);
    plot(t,Isyntrain-3,'-','Color',[.7 .7 .7],'linewidth',2);
    plot(tsepeak,vsemax,'or')
    plot(tsetrough,vsemin,'or')
    axis([0 Tmax -4 5]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('');
    legend('Current-based');
    
%     figure
%     hold on
%     plot(CycleFreqOrdPeak,vsemaxord,'ob','linewidth',2);
%     plot(CycleFreqOrdTrough,vseminord,'ob','linewidth',2);
%     axis([0 100 -2 4]);
%     set(gca,'fontsize',24);
%     xlabel('f  [Hz]');
%     ylabel('v_{max}  [mV]');
    
    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(CycleFreqOrd2,vsemaxord2,'ob','linewidth',2);
    plot(CycleFreqOrd2,vseminord2,'ob','linewidth',2);
    axis([0 100 -2 4]);
    set(gca,'fontsize',24);
    xlabel('f  [Hz]');
    ylabel('v_{max}  [mV]');
     legend('Current-based');

%     % Inhibitory synaptic-like current-based inputs
% 
%     Gsyn = -1;
% 
%     % Computation of the numerical solution
% 
%     vs = zeros(1,length(t));
%     ws = zeros(1,length(t));
% 
%     for j=1:length(t)-1
%         k1v = -gL*vs(j)-g1*ws(j)+Gsyn*Isyntrain(j);
%         k1w = (vs(j)-ws(j))/tau1;
%         av = vs(j)+k1v*dt;
%        `av = av+sqrt(2*D*dt)*eta(j);
%         aw = ws(j)+k1w*dt;
%         k2v = -gL*av-g1*aw+Gsyn*Isyntrain(j+1);
%         k2w = (av-aw)/tau1;
%         vs(j+1) = vs(j)+(k1v+k2v)*dt/2;
%        vs(j+1) = vs(j+1)+sqrt(2*D*dt)*eta(j);
%         ws(j+1) = ws(j)+(k1w+k2w)*dt/2;
%     end
% 
% 
%     % Computation of the peaks and peak-times
% 
%     [tpeakzap,zapmax,contpeakzap] = PeaksOsc(Isyntrain,t,dt,Tmax,dt);
%     tpeakzap = [0 tpeakzap];
%     zapmax = [1 zapmax];
%     [tsipeak,vsimax,cntpeak] = PeaksOsc(vs,t,dt,Tmax,dt);
%     [tsitrough,vsimin,cnttrough] = TroughsOsc(vs,t,dt,Tmax,dt);
%     tsipeak = [tsipeak t(end)];
%     vsimax = [vsimax vs(end)];
%     Phases = zeros(1,Ncycles);
%     for k=1:Ncycles
%         Phases(k) = (tsipeak(k)-tpeakzap(k))*CycleFreq(k);
%     end
%     vsimax2 = zeros(1,Ncycles);
%     vsimin2 = zeros(1,Ncycles);
%     tsipeak2 = zeros(1,Ncycles);
%     tsitrough2 = zeros(1,Ncycles);
%     for k=1:Ncycles
%         [vsimax2(k),jpeak2(k)] = max(vs(jbdry(k):jbdry(k+1)-1));  
%         [vsimin2(k),jtrough2(k)] = min(vs(jbdry(k):jbdry(k+1)-1)); 
%         tsipeak2(k) = tbdry(k)+jpeak2(k)*dt;
%         tsitrough2(k) = tbdry(k)+jtrough2(k)*dt;
%     end
%     [CycleFreqOrd2,vsimaxord2,vsiminord2] = OrderingCycleFreqPeakTrough(CycleFreq,vsimax2,vsimin2);
%     CycleFreqOutputPeak = 1000./diff(tsitrough);
%     [CycleFreqOrdPeak,vsimaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsimax(1:end-1));
%     CycleFreqOutputTrough = 1000./diff(tsipeak);
%     [CycleFreqOrdTrough,vsiminord] = OrderingCycleFreqTrough(CycleFreqOutputTrough,vsimin(2:end));     
%     
% 
%     vsi = vs;
    

    % figure
    % hold on
    % plot(t,Isyntrain,'-b','linewidth',2);
    % axis([0 Tmax -1.2 1.2]);
    % set(gca,'fontsize',24);
    % xlabel('t');
    % ylabel('');

%     figure
%     hold on
%     plot(t,vsi,'-b','linewidth',2);
%     plot(t,Isyntrain-5,'-','Color',[.7 .7 .7],'linewidth',2);
%     plot(tsipeak,vsimax,'or')
%     plot(tsitrough,vsimin,'or')
%     axis([0 Tmax -6 4]);
%     set(gca,'fontsize',24);
%     xlabel('t  [ms]');
%     ylabel('');
%     
% %     figure
% %     hold on
% %     plot(CycleFreqOrdPeak,vsimaxord,'ob','linewidth',2);
% %     plot(CycleFreqOrdTrough,vsiminord,'ob','linewidth',2);
% %     axis([0 100 -2 4]);
% %     set(gca,'fontsize',24);
% %     xlabel('f  [Hz]');
% %     ylabel('v_{max}  [mV]');
%     
%     figure
%     hold on
%     plot(-100,-100,'ob','linewidth',2);
%     plot(-100,-100,'or','linewidth',2);
%     plot(CycleFreqOrd2,vsimaxord2,'ob','linewidth',2);
%     plot(CycleFreqOrd2,vsiminord2,'ob','linewidth',2);
%     axis([0 100 -4 2]);
%     set(gca,'fontsize',24);
%     xlabel('f  [Hz]');
%     ylabel('v_{max}  [mV]');

    




    % % Conductance-based synaptic-like inputs
    % 
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
    
    % Computation of the peaks and peak-times
    
    vscbmax = zeros(1,Ncycles);
    jpeak = zeros(1,Ncycles);
    tpeak = zeros(1,Ncycles);
    Phasescb = zeros(1,Ncycles);
    vscbmin = zeros(1,Ncycles);
    ttrough = zeros(1,Ncycles);
    jtrough = zeros(1,Ncycles);
    for k=1:Ncycles
        [vscbmax(k),jpeak(k)] = max(vscb(jbdry(k):jbdry(k+1)));    
        tpeak(k) = tbdry(k)+jpeak(k)*dt;
        Phasescb(k) = (tpeak(k)-tbdry(k)-CyclePer(k)/2)*CycleFreq(k);   
        [vscbmin(k),jtrough(k)] = min(vscb(jbdry(k):jbdry(k+1))); 
        ttrough(k) = tbdry(k)+jtrough(k)*dt;
    end
    Phasescb = Phasescb*2*pi/1000;
    
    
    [CycleFreqOrd,vscbmaxord,phasescbord,vscbminord] = OrderingCycleFreq(CycleFreq,vscbmax,Phasescb,vscbmin);
    
    vscbe = vscb;
    vscbemaxord = vscbmaxord;
    vscbeminord = vscbminord;
    tscbepeak = tpeak;
    tscbetrough = ttrough;
    
    % figure
    % hold on
    % plot(t,Gex*Isyntrain.*(vscb-Ex),'-b','linewidth',2);
    % axis([0 Tmax -1.2 1.2]);
    % set(gca,'fontsize',24);
    % xlabel('t');
    % ylabel('');
    
    Izaplikecbe = -(Gsyncb*Isyntrain.*(vscbe-Esyn));
    figure
    hold on
    plot(t,vscbe,'-b','linewidth',2);
    %plot(t,-Isyntrain.*(vscb-Eex)/max(-Isyntrain.*(vscb-Eex))-2,'-','Color',[.7 .7 .7],'linewidth',2);
    plot(tscbepeak,vscbmax,'or','linewidth',2)
    plot(tscbetrough,vscbmin,'or','linewidth',2)
    %plot(t,Izaplikecbe/max(Izaplikecbe)-2,'-','Color',[.7 .7 .7],'linewidth',2);
    plot(t,Izaplikecbe-2,'-','Color',[.7 .7 .7],'linewidth',2);
    axis([0 Tmax -3 5]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('');
     legend('Conductance-based');
    
    figure
    hold on
    plot(CycleFreqOrd,vscbemaxord,'ob','linewidth',2);
    plot(CycleFreqOrd,vscbeminord,'ob','linewidth',2);
    %plot(CycleFreqOrd,(vscbemaxord-vscbeminord)/2,'or','linewidth',2);
    axis([0 100 -3 2]);
    set(gca,'fontsize',24);
    xlabel('f  [Hz]');
    ylabel('v_{max}  [mV]');
    legend('Conductance-based');
    
    % % Inhibitory
    % 
    % Esyn = -1;
    % 
    % % Computation of the numerical solution
    % 
    % vscb = zeros(1,length(t));
    % wscb = zeros(1,length(t));
    % 
    % for j=1:length(t)-1
    %     k1v = -gL*vscb(j)-g1*wscb(j)-Gsyncb*Isyntrain(j)*(vscb(j)-Esyn);
    %     k1w = (vscb(j)-wscb(j))/tau1;
    %     av = vscb(j)+k1v*dt;
    %     av = av+sqrt(2*D*dt)*eta(j);
    %     aw = wscb(j)+k1w*dt;
    %     k2v = -gL*av-g1*aw-Gsyncb*Isyntrain(j+1)*(av-Esyn);
    %     k2w = (av-aw)/tau1;
    %     vscb(j+1) = vscb(j)+(k1v+k2v)*dt/2;
    %     vscb(j+1) = vs(j+1)+sqrt(2*D*dt)*eta(j);
    %     wscb(j+1) = wscb(j)+(k1w+k2w)*dt/2;
    % end
    % 
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
    % 
    % [CycleFreqOrd,vscbmaxord,phasescbord,vscbminord] = OrderingCycleFreq(CycleFreq,vscbmax,Phasescb,vscbmin);
    % 
    % vscbi = vscb;
    % vscbimaxord = vscbmaxord;
    % vscbiminord = vscbminord;
    % tscbipeak = tpeak;
    % tscbitrough = ttrough;
    % 
    % % figure
    % % hold on
    % % plot(t,Gex*Isyntrain.*(vscb-Ex),'-b','linewidth',2);
    % % axis([0 Tmax -1.2 1.2]);
    % % set(gca,'fontsize',24);
    % % xlabel('t');
    % % ylabel('');
    % 
    % figure
    % hold on
    % plot(t,vscbi,'-b','linewidth',2);
    % %plot(t,-Isyntrain.*(vscb-Eex)/max(-Isyntrain.*(vscb-Eex))-2,'-','Color',[.7 .7 .7],'linewidth',2);
    % plot(t,(Gsyncb*Isyntrain.*(vscbi-Esyn))/(max(abs(-Gsyncb*Isyntrain.*(vscbi-Esyn))))-2,'-','Color',[.7 .7 .7],'linewidth',2);
    % axis([0 Tmax -3 5]);
    % set(gca,'fontsize',24);
    % xlabel('t  [ms]');
    % ylabel('');
    % 
    % figure
    % hold on
    % plot(CycleFreqOrd,vscbimaxord,'ob','linewidth',2);
    % plot(CycleFreqOrd,vscbiminord,'ob','linewidth',2);
    % plot(CycleFreqOrd,(vscbimaxord-vscbiminord)/2,'or','linewidth',2);
    % axis([0 100 -3 2]);
    % set(gca,'fontsize',24);
    % xlabel('f  [Hz]');
    % ylabel('v_{max}  [mV]');
    % 
    
elseif SYN == 2
    
    Ntrials = 100;
    
    % Sequence of synaptic-like inputs
    % First one: regular order
    % Others: random order

    % Generation of the ZAP-like regular and random synaptic-like funtions

    taurse_e = 1;     % AMPA
    taudec_e = 5;     % AMPA
    taurse_i = 0.2;   % GABA_A     
    taudec_i = 10;    % GABA_A
    Eex = 20;
    Ein = -60;

    taudec = taudec_e;

    
    % Generation of the ZAP-like regular and random cycle periods

    RNDM = 0;
    if RNDM == 0
        CycleFreq = 1:1:Ncycles;
    elseif RNDM == 1
        CycleFreq = randperm(Ncycles,Ncycles);
    end
    CyclePer = 1000./CycleFreq;


    t = 0:dt:CyclePer(1);
    Isyntrain = Fsyn(t,Ain,t(1),taudec);
    tbdry(1) = 0;
    jbdry(1) = 1;
    tbdry(2) = t(end);
    jbdry(2) = length(t);
    for k=2:Ncycles
        tloc = dt:dt:CyclePer(k);  
        Floc = Fsyn(tloc,Ain,tloc(1),taudec);
        t = [t t(end)+tloc];
        Isyntrain = [Isyntrain Floc];
        tbdry(k+1) = t(end);
        jbdry(k+1) = length(t);
    end
    Tmax = t(end);
    Deltajbdry = diff(jbdry);
    
    % White (Gaussian) noise

    eta = randn(1,length(t));

    % Excitatory synaptic-like current-based inputs

    Gsyn = 1;

    % Computation of the numerical solution

    vs = zeros(1,length(t));
    ws = zeros(1,length(t));

    for j=1:length(t)-1
        k1v = -gL*vs(j)-g1*ws(j)+Gsyn*Isyntrain(j);
        k1w = (vs(j)-ws(j))/tau1;
        av = vs(j)+k1v*dt;
        av = av+sqrt(2*D*dt)*eta(j);
        aw = ws(j)+k1w*dt;
        k2v = -gL*av-g1*aw+Gsyn*Isyntrain(j+1);
        k2w = (av-aw)/tau1;
        vs(j+1) = vs(j)+(k1v+k2v)*dt/2;
        vs(j+1) = vs(j+1)+sqrt(2*D*dt)*eta(j);
        ws(j+1) = ws(j)+(k1w+k2w)*dt/2;
    end


    % Computation of the peaks and peak-times

    [tpeakzap,zapmax,contpeakzap] = PeaksOsc(Isyntrain,t,dt,Tmax,dt);
    tpeakzap = [0 tpeakzap];
    zapmax = [1 zapmax];
    [tsepeak,vsemax,cntpeak] = PeaksOsc(vs,t,dt,Tmax,dt);
    [tsetrough,vsemin,cnttrough] = TroughsOsc(vs,t,dt,Tmax,dt);
    tsetrough = [tsetrough t(end)];
    vsemin = [vsemin vs(end)];
    Phases = zeros(1,Ncycles);
    for k=1:Ncycles
        Phases(k) = (tsepeak(k)-tpeakzap(k))*CycleFreq(k);
    end
    vsemax2 = zeros(1,Ncycles);
    vsemin2 = zeros(1,Ncycles);
    tsepeak2 = zeros(1,Ncycles);
    tsetrough2 = zeros(1,Ncycles);
    for k=1:Ncycles
        [vsemax2(k),jpeak2(k)] = max(vs(jbdry(k):jbdry(k+1)-1));  
        [vsemin2(k),jtrough2(k)] = min(vs(jbdry(k):jbdry(k+1)-1)); 
        tsepeak2(k) = tbdry(k)+jpeak2(k)*dt;
        tsetrough2(k) = tbdry(k)+jtrough2(k)*dt;
    end
    [CycleFreqOrd2,vsemaxord2,vseminord2] = OrderingCycleFreqPeakTrough(CycleFreq,vsemax2,vsemin2);
    CycleFreqOutputPeak = 1000./diff(tsetrough);
    if length(CycleFreqOutputPeak) < length(vsemax)   
        [CycleFreqOrdPeak,vsemaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsemax(1:end-1));
    else
        [CycleFreqOrdPeak,vsemaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsemax);
    end
    CycleFreqOutputTrough = 1000./diff(tsepeak);
   [CycleFreqOrdTrough,vseminord] = OrderingCycleFreqTrough(CycleFreqOutputTrough,vsemin(2:end));     
 

    vse = vs;
    
    vse_scale = vsemax(1)-vsemin(1);
    
    figure
    hold on
    plot(t,vse,'-b','linewidth',2);
    plot(t,Isyntrain-3,'-','Color',[.7 .7 .7],'linewidth',2);
    plot(tsepeak,vsemax,'or')
    plot(tsetrough,vsemin,'or')
    axis([0 Tmax -4 5]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('');
    legend('Current-based');
    
%     figure
%     hold on
%     plot(-100,-100,'ob','linewidth',2);
%     plot(-100,-100,'or','linewidth',2);
%     plot(CycleFreqOrd2,vsemaxord2,'ob','linewidth',2);
%     plot(CycleFreqOrd2,vseminord2,'ob','linewidth',2);
%     axis([0 100 -2 4]);
%     set(gca,'fontsize',24);
%     xlabel('f  [Hz]');
%     ylabel('v_{max}  [mV]');
%     title('Current-based');
    
     % Conductance-based synaptic-like inputs
     
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
    
    % Computation of the peaks and peak-times
    
    vscbmax = zeros(1,Ncycles);
    jpeak = zeros(1,Ncycles);
    tpeak = zeros(1,Ncycles);
    Phasescb = zeros(1,Ncycles);
    vscbmin = zeros(1,Ncycles);
    ttrough = zeros(1,Ncycles);
    jtrough = zeros(1,Ncycles);
    for k=1:Ncycles
        [vscbmax(k),jpeak(k)] = max(vscb(jbdry(k):jbdry(k+1)));    
        tpeak(k) = tbdry(k)+jpeak(k)*dt;
        Phasescb(k) = (tpeak(k)-tbdry(k)-CyclePer(k)/2)*CycleFreq(k);   
        [vscbmin(k),jtrough(k)] = min(vscb(jbdry(k):jbdry(k+1))); 
        ttrough(k) = tbdry(k)+jtrough(k)*dt;
    end
    Phasescb = Phasescb*2*pi/1000;
    
    
    [CycleFreqOrd,vscbmaxord,phasescbord,vscbminord] = OrderingCycleFreq(CycleFreq,vscbmax,Phasescb,vscbmin);
    
    vscbe = vscb;
    vscbemaxord = vscbmaxord;
    vscbeminord = vscbminord;
    tscbepeak = tpeak;
    tscbetrough = ttrough;
    
    Ampse_bse = vsemax(1)-vsemin(1);
    Ampscbe_bse = vscbmax(1)-vscbmin(1);    
    
    
    Izaplikecbe = -(Gsyncb*Isyntrain.*(vscbe-Esyn));
    figure
    hold on
    plot(t,vscbe,'-b','linewidth',2);
    %plot(t,-Isyntrain.*(vscb-Eex)/max(-Isyntrain.*(vscb-Eex))-2,'-','Color',[.7 .7 .7],'linewidth',2);
    plot(tscbepeak,vscbmax,'or','linewidth',2)
    plot(tscbetrough,vscbmin,'or','linewidth',2)
    %plot(t,Izaplikecbe/max(Izaplikecbe)-2,'-','Color',[.7 .7 .7],'linewidth',2);
    plot(t,Izaplikecbe-2,'-','Color',[.7 .7 .7],'linewidth',2);
    axis([0 Tmax -3 5]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('');
     legend('Conductance-based');
    
%     figure
%     hold on
%     plot(CycleFreqOrd,vscbemaxord,'ob','linewidth',2);
%     plot(CycleFreqOrd,vscbeminord,'ob','linewidth',2);
%     %plot(CycleFreqOrd,(vscbemaxord-vscbeminord)/2,'or','linewidth',2);
%     axis([0 100 -3 2]);
%     set(gca,'fontsize',24);
%     xlabel('f  [Hz]');
%     ylabel('v_{max}  [mV]');
%     title('Conductance-based');
    
    
    vsemax_bse = vsemaxord2;
    vsemin_bse = vseminord2;
    vsemaxord_vec = zeros(Ntrials,Ncycles);
    vseminord_vec = zeros(Ntrials,Ncycles);
    vscbemax_bse = vscbemaxord;
    vscbemin_bse = vscbeminord;
    vscbemaxord_vec = zeros(Ntrials,Ncycles);
    vscbeminord_vec = zeros(Ntrials,Ncycles);
    
    vs_vec = zeros(Ntrials,length(t));
    vscb_vec = zeros(Ntrials,length(t));
    
    for l=1:Ntrials
       
        l
       
        RNDM = 1;
        if RNDM == 0
            CycleFreq = 1:1:Ncycles;
        elseif RNDM == 1
            CycleFreq = randperm(Ncycles,Ncycles);
        end
        CyclePer = 1000./CycleFreq;
       
        t = 0:dt:CyclePer(1);
        Isyntrain = Fsyn(t,Ain,t(1),taudec);
        tbdry(1) = 0;
        jbdry(1) = 1;
        tbdry(2) = t(end);
        jbdry(2) = length(t);
        for k=2:Ncycles
            tloc = dt:dt:CyclePer(k);  
            Floc = Fsyn(tloc,Ain,tloc(1),taudec);
            t = [t t(end)+tloc];
            Isyntrain = [Isyntrain Floc];
            tbdry(k+1) = t(end);
            jbdry(k+1) = length(t);
        end
        Tmax = t(end);
        Deltajbdry = diff(jbdry);
        
        % Current-based inputs
        
        % Computation of the numerical solution

        vs = zeros(1,length(t));
        ws = zeros(1,length(t));

        for j=1:length(t)-1
            k1v = -gL*vs(j)-g1*ws(j)+Gsyn*Isyntrain(j);
            k1w = (vs(j)-ws(j))/tau1;
            av = vs(j)+k1v*dt;
            av = av+sqrt(2*D*dt)*eta(j);
            aw = ws(j)+k1w*dt;
            k2v = -gL*av-g1*aw+Gsyn*Isyntrain(j+1);
            k2w = (av-aw)/tau1;
            vs(j+1) = vs(j)+(k1v+k2v)*dt/2;
            vs(j+1) = vs(j+1)+sqrt(2*D*dt)*eta(j);
            ws(j+1) = ws(j)+(k1w+k2w)*dt/2;
        end

        vs_vec(l,:) = vs;

        % Computation of the peaks and peak-times

        [tpeakzap,zapmax,contpeakzap] = PeaksOsc(Isyntrain,t,dt,Tmax,dt);
        tpeakzap = [0 tpeakzap];
        zapmax = [1 zapmax];
        [tsepeak,vsemax,cntpeak] = PeaksOsc(vs,t,dt,Tmax,dt);
        [tsetrough,vsemin,cnttrough] = TroughsOsc(vs,t,dt,Tmax,dt);
        tsetrough = [tsetrough t(end)];
        vsemin = [vsemin vs(end)];
        Phases = zeros(1,Ncycles);
        for k=1:Ncycles
            Phases(k) = (tsepeak(k)-tpeakzap(k))*CycleFreq(k);
        end
        vsemax2 = zeros(1,Ncycles);
        vsemin2 = zeros(1,Ncycles);
        tsepeak2 = zeros(1,Ncycles);
        tsetrough2 = zeros(1,Ncycles);
        for k=1:Ncycles
            [vsemax2(k),jpeak2(k)] = max(vs(jbdry(k):jbdry(k+1)-1));  
            [vsemin2(k),jtrough2(k)] = min(vs(jbdry(k):jbdry(k+1)-1)); 
            tsepeak2(k) = tbdry(k)+jpeak2(k)*dt;
            tsetrough2(k) = tbdry(k)+jtrough2(k)*dt;
        end
        [CycleFreqOrd2,vsemaxord2,vseminord2] = OrderingCycleFreqPeakTrough(CycleFreq,vsemax2,vsemin2);
        CycleFreqOutputPeak = 1000./diff(tsetrough);
        if length(CycleFreqOutputPeak) < length(vsemax)   
            [CycleFreqOrdPeak,vsemaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsemax(1:end-1));
        else
            [CycleFreqOrdPeak,vsemaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsemax);
        end
        CycleFreqOutputTrough = 1000./diff(tsepeak);
        [CycleFreqOrdTrough,vseminord] = OrderingCycleFreqTrough(CycleFreqOutputTrough,vsemin(2:end));     

        vsemaxord_vec(l,1:length(vsemaxord2)) = vsemaxord2;
        vseminord_vec(l,1:length(vseminord2)) = vseminord2;
    
        % Conductance-based synaptic-like inputs
        
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
        
        vscb_vec(l,:) = vscb;

        % Computation of the peaks and peak-times

        vscbmax = zeros(1,Ncycles);
        jpeak = zeros(1,Ncycles);
        tpeak = zeros(1,Ncycles);
        Phasescb = zeros(1,Ncycles);
        vscbmin = zeros(1,Ncycles);
        ttrough = zeros(1,Ncycles);
        jtrough = zeros(1,Ncycles);
        for k=1:Ncycles
            [vscbmax(k),jpeak(k)] = max(vscb(jbdry(k):jbdry(k+1)));    
            tpeak(k) = tbdry(k)+jpeak(k)*dt;
            Phasescb(k) = (tpeak(k)-tbdry(k)-CyclePer(k)/2)*CycleFreq(k);   
            [vscbmin(k),jtrough(k)] = min(vscb(jbdry(k):jbdry(k+1))); 
            ttrough(k) = tbdry(k)+jtrough(k)*dt;
        end
        Phasescb = Phasescb*2*pi/1000;


        [CycleFreqOrd,vscbmaxord,phasescbord,vscbminord] = OrderingCycleFreq(CycleFreq,vscbmax,Phasescb,vscbmin);

        vscbe = vscb;
        vscbemaxord = vscbmaxord;
        vscbeminord = vscbminord;
        tscbepeak = tpeak;
        tscbetrough = ttrough;
    
        vscbemaxord_vec(l,1:length(vscbemaxord)) = vscbemaxord;
        vscbeminord_vec(l,1:length(vscbemaxord)) = vscbeminord;
       
    end
    
    
   figure
   hold on
   for l=1:Ntrials
       plot(vsemaxord_vec(l,:),'o')
       plot(vseminord_vec(l,:),'o')
   end
   plot(vsemax_bse,'ob','linewidth',2)
   plot(vsemin_bse,'ob','linewidth',2)
   axis([0 100 -2 4]);
   set(gca,'fontsize',24);
   xlabel('f  [Hz]');
   ylabel('V_{ENV}^{+/-}  [mV]');
   title('Current-based');
    
   figure
   hold on
   for l=1:Ntrials
       plot(vscbemaxord_vec(l,:),'o')
       plot(vscbeminord_vec(l,:),'o')
   end
   plot(vscbemax_bse,'ob','linewidth',2)
   plot(vscbemin_bse,'ob','linewidth',2)
   axis([0 100 -2 4]);
   set(gca,'fontsize',24);
   xlabel('f  [Hz]');
   ylabel('V_{ENV}^{+/-}  [mV]');
   title('Conductance-based');
   

   vsemax_mean = zeros(1,Ncycles);
   vsemin_mean = zeros(1,Ncycles);
   vsemax_var = zeros(1,Ncycles);
   vsemin_var = zeros(1,Ncycles);
   for k=1:Ncycles
       vsemax_mean(k) = mean(vsemaxord_vec(:,k));
       vsemin_mean(k) = mean(vseminord_vec(:,k));
       vsemax_var(k) = var(vsemaxord_vec(:,k));
       vsemin_var(k) = var(vseminord_vec(:,k));
   end
        
    vscbemax_mean = zeros(1,Ncycles);
    vscbemin_mean = zeros(1,Ncycles);
    vscbemax_var = zeros(1,Ncycles);
    vscbemin_var = zeros(1,Ncycles);
    for k=1:Ncycles
         vscbemax_mean(k) = mean(vscbemaxord_vec(:,k));
         vscbemin_mean(k) = mean(vscbeminord_vec(:,k));
         vscbemax_var(k) = var(vscbemaxord_vec(:,k));
         vscbemin_var(k) = var(vscbeminord_vec(:,k));
    end
      
   figure
   hold on
   plot(vsemax_var,'ob','linewidth',2) 
   plot(vsemin_var,'or','linewidth',2)
   plot(vsemax_var./Ampse_bse,'o','Color',lightblueish,'linewidth',2); 
   plot(vsemin_var./Ampse_bse,'o','Color',lightcoral,'linewidth',2);
   set(gca,'fontsize',24);
   xlabel('f  [Hz]');
   ylabel('Var');
   title('Current-based');   
   legend('V_{ENV}^+','V_{ENV}^-','V_{ENV,norm}^+','V_{ENV,norm}^-','Location','NorthWest');
   
   figure
   hold on
   plot(vscbemax_var,'ob','linewidth',2) 
   plot(vscbemin_var,'or','linewidth',2)
   plot(vscbemax_var./Ampscbe_bse,'o','Color',lightblueish,'linewidth',2) 
   plot(vscbemin_var./Ampscbe_bse,'o','Color',lightcoral,'linewidth',2)
   set(gca,'fontsize',24);
   xlabel('f  [Hz]');
   ylabel('Var (V_{ENV}^{\pm})');
   title('Conductance-based');
   legend('V_{ENV}^+','V_{ENV}^-','V_{ENV,norm}^+','V_{ENV,norm}^-','Location','NorthWest');
   
   figure
   hold on
   plot(vsemax_var,'ob','linewidth',2) 
   plot(vsemin_var,'or','linewidth',2)
   plot(vscbemax_var,'o','Color',lightblueish,'linewidth',2); 
   plot(vscbemin_var,'o','Color',lightcoral,'linewidth',2);
   set(gca,'fontsize',24);
   xlabel('f  [Hz]');
   ylabel('Var (V_{ENV}^{\pm})'); 
   legend('V_{ENV,curr}^+','V_{ENV,curr}^-' ,'V_{ENV,cond}^+','V_{ENV,cond}^-','Location','NorthWest');
   
   figure
   hold on
   plot(vsemax_var./Ampse_bse,'ob','linewidth',2) 
   plot(vsemin_var./Ampse_bse,'or','linewidth',2)
   plot(vscbemax_var./Ampscbe_bse,'o','Color',lightblueish,'linewidth',2); 
   plot(vscbemin_var./Ampscbe_bse,'o','Color',lightcoral,'linewidth',2);
   set(gca,'fontsize',24);
   xlabel('f  [Hz]');
   ylabel('              Var (V_{ENV}^{\pm})    (norm)'); 
   legend('V_{ENV,curr}^+','V_{ENV,curr}^-' ,'V_{ENV,cond}^+','V_{ENV,cond}^-','Location','NorthWest');
  
   
%    Variab_vsemax = zeros(1,Ncycles);
%    Variab_vsemin = zeros(1,Ncycles);
%    for k=1:Ncycles
%        Variab_vsemax(k) = sum((vsemaxord_vec(:,k)-vsemax_bse(k)).^2);
%        Variab_vsemin(k) = sum((vseminord_vec(:,k)-vsemin_bse(k)).^2);
%    end

%     Variab_vscbemax = zeros(1,Ncycles);
%     Variab_vscbemin = zeros(1,Ncycles);
%     for k=1:Ncycles
%        Variab_vscbemax(k) = sum((vscbemaxord_vec(:,k)-vscbemax_bse(k)).^2);
%        Variab_vscbemin(k) = sum((vscbeminord_vec(:,k)-vscbemin_bse(k)).^2);
%     end
        
elseif SYN == 3
    
    % Not completed since it does not seem to be useful
    
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

    
    % Generation of the ZAP-like regular and random synaptic-like funtions

    taurse_e = 1;     % AMPA
    taudec_e = 5;     % AMPA
    taurse_i = 0.2;   % GABA_A     
    taudec_i = 10;    % GABA_A
    Eex = 20;
    Ein = -60;

    taudec = taudec_e;

    
    % Generation of the ZAP-like regular and random cycle periods

    RNDM = 0;
    if RNDM == 0
        CycleFreq = 1:1:Ncycles;
    elseif RNDM == 1
        CycleFreq = randperm(Ncycles,Ncycles);
    end
    CyclePer = 1000./CycleFreq;
    
    
    taudec = taudec_e;
    

    t = 0:dt:CyclePer(1);
    Isyntrain = Fsyn(t,Ain,t(1),taudec);
    tbdry(1) = 0;
    jbdry(1) = 1;
    tbdry(2) = t(end);
    jbdry(2) = length(t);
    for k=2:Ncycles
        tloc = dt:dt:CyclePer(k);  
        Floc = Fsyn(tloc,Ain,tloc(1),taudec);
        t = [t t(end)+tloc];
        Isyntrain = [Isyntrain Floc];
        tbdry(k+1) = t(end);
        jbdry(k+1) = length(t);
    end
    Tmax = t(end);
    Deltajbdry = diff(jbdry);
    
    % White (Gaussian) noise

    eta = randn(1,length(t));
    
    % Excitatory synaptic-like current-based inputs

    Gsyn = 2;

    % Computation of the numerical solution: linear system

    vs = zeros(1,length(t));
    ws = zeros(1,length(t));

    for j=1:length(t)-1
        k1v = -gL*vs(j)-g1*ws(j)+Gsyn*Isyntrain(j);
        k1w = (vs(j)-ws(j))/tau1;
        av = vs(j)+k1v*dt;
        av = av+sqrt(2*D*dt)*eta(j);
        aw = ws(j)+k1w*dt;
        k2v = -gL*av-g1*aw+Gsyn*Isyntrain(j+1);
        k2w = (av-aw)/tau1;
        vs(j+1) = vs(j)+(k1v+k2v)*dt/2;
        vs(j+1) = vs(j+1)+sqrt(2*D*dt)*eta(j);
        ws(j+1) = ws(j)+(k1w+k2w)*dt/2;
    end


    % Computation of the peaks and peak-times

    [tpeakzap,zapmax,contpeakzap] = PeaksOsc(Isyntrain,t,dt,Tmax,dt);
    tpeakzap = [0 tpeakzap];
    zapmax = [1 zapmax];
    [tsepeak,vsemax,cntpeak] = PeaksOsc(vs,t,dt,Tmax,dt);
    [tsetrough,vsemin,cnttrough] = TroughsOsc(vs,t,dt,Tmax,dt);
    tsetrough = [tsetrough t(end)];
    vsemin = [vsemin vs(end)];
    Phases = zeros(1,Ncycles);
    for k=1:Ncycles
        Phases(k) = (tsepeak(k)-tpeakzap(k))*CycleFreq(k);
    end
    vsemax2 = zeros(1,Ncycles);
    vsemin2 = zeros(1,Ncycles);
    tsepeak2 = zeros(1,Ncycles);
    tsetrough2 = zeros(1,Ncycles);
    for k=1:Ncycles
        [vsemax2(k),jpeak2(k)] = max(vs(jbdry(k):jbdry(k+1)-1));  
        [vsemin2(k),jtrough2(k)] = min(vs(jbdry(k):jbdry(k+1)-1)); 
        tsepeak2(k) = tbdry(k)+jpeak2(k)*dt;
        tsetrough2(k) = tbdry(k)+jtrough2(k)*dt;
    end
    [CycleFreqOrd2,vsemaxord2,vseminord2] = OrderingCycleFreqPeakTrough(CycleFreq,vsemax2,vsemin2);
    CycleFreqOutputPeak = 1000./diff(tsetrough);
    if length(CycleFreqOutputPeak) < length(vsemax)   
        [CycleFreqOrdPeak,vsemaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsemax(1:end-1));
    else
        [CycleFreqOrdPeak,vsemaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsemax);
    end
    CycleFreqOutputTrough = 1000./diff(tsepeak);
    [CycleFreqOrdTrough,vseminord] = OrderingCycleFreqTrough(CycleFreqOutputTrough,vsemin(2:end));     
 
    tsepeaklin = tsepeak;
    tsetroughlin = tsetrough;
    vsemaxlin = vsemax;
    vseminlin = vsemin;
    vselin = vs;
    wselin = ws;
    
    
  
    
    % % Computation of the numerical solution: PWL system

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


    % Computation of the peaks and peak-times

    [tpeakzap,zapmax,contpeakzap] = PeaksOsc(Isyntrain,t,dt,Tmax,dt);
    tpeakzap = [0 tpeakzap];
    zapmax = [1 zapmax];
    [tsepeak,vsemax,cntpeak] = PeaksOsc(vs,t,dt,Tmax,dt);
    [tsetrough,vsemin,cnttrough] = TroughsOsc(vs,t,dt,Tmax,dt);
    tsetrough = [tsetrough t(end)];
    vsemin = [vsemin vs(end)];
    Phases = zeros(1,Ncycles);
    for k=1:Ncycles
        Phases(k) = (tsepeak(k)-tpeakzap(k))*CycleFreq(k);
    end
    vsemax2 = zeros(1,Ncycles);
    vsemin2 = zeros(1,Ncycles);
    tsepeak2 = zeros(1,Ncycles);
    tsetrough2 = zeros(1,Ncycles);
    for k=1:Ncycles
        [vsemax2(k),jpeak2(k)] = max(vs(jbdry(k):jbdry(k+1)-1));  
        [vsemin2(k),jtrough2(k)] = min(vs(jbdry(k):jbdry(k+1)-1)); 
        tsepeak2(k) = tbdry(k)+jpeak2(k)*dt;
        tsetrough2(k) = tbdry(k)+jtrough2(k)*dt;
    end
    [CycleFreqOrd2,vsemaxord2,vseminord2] = OrderingCycleFreqPeakTrough(CycleFreq,vsemax2,vsemin2);
    CycleFreqOutputPeak = 1000./diff(tsetrough);
    if length(CycleFreqOutputPeak) < length(vsemax)   
        [CycleFreqOrdPeak,vsemaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsemax(1:end-1));
    else
        [CycleFreqOrdPeak,vsemaxord] = OrderingCycleFreqPeak(CycleFreqOutputPeak,vsemax);
    end
    CycleFreqOutputTrough = 1000./diff(tsepeak);
    [CycleFreqOrdTrough,vseminord] = OrderingCycleFreqTrough(CycleFreqOutputTrough,vsemin(2:end));     
 
    tsepeakpwl = tsepeak;
    tsetroughpwl = tsetrough;
    vsemaxpwl = vsemax;
    vseminpwl = vsemin;
    vsepwl = vs;
    wsepwl = ws;
    
    figure
    hold on
    plot(t,vsepwl,'-b','linewidth',2);
    plot(t,Isyntrain-3,'-','Color',[.7 .7 .7],'linewidth',2);
    plot(tsepeakpwl,vsemaxpwl,'or')
    plot(tsetroughpwl,vseminpwl,'or')
    plot(tsepeakpwl,vsemaxlin,'o','Color',lightcoral,'linewidth',2)
    plot(tsetroughpwl,vseminlin,'o','Color',lightcoral,'linewidth',2)
    %plot(t,vselin,'--g','linewidth',2);
    axis([0 Tmax -4 5]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('');
    legend('Current-based');
    
    figure(101)
    plot(vsepwl,wselin,'-','Color',lightblueish,'linewidth',1);
    plot(vselin,wselin,'--','Color',lightcoral,'linewidth',1);
    legend('v-nullclne','w-nullcline','trajectory','Location','NorthEast');

end

PWS = 3;
if PWS == 1

    % Power spectrum density computation
    
    if SIN == 1

   
        
        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,Izaplike);
        ZAP_Psd = PsdManSmooth;

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,v);
        V_Psd = PsdManSmooth;

        Z = V_Psd./ZAP_Psd;

        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,Z,'ob','linewidth',1);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,20]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('Z');
        
        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,V_Psd,'ob','linewidth',1);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,2]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('PSD');
        
    elseif SQW == 1
        
        
        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,Isqwtrain);
        ZAP_Psd = PsdManSmooth;

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vsqw);
        V_Psd = PsdManSmooth;

        Z = V_Psd./ZAP_Psd;

        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,Z,'ob','linewidth',1);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,20]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('Z');
        
         figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,V_Psd,'ob','linewidth',1);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,2]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('PSD');
        
    elseif SYN == 1
        

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,Isyntrain);
        ZAP_Psd = PsdManSmooth;

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vse);
        V_Psd = PsdManSmooth;

        Ze = V_Psd./ZAP_Psd;

        figure
        hold on
        plot(freqbin,Ze,'ob','linewidth',1);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,20]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('Z');
        legend('Current-based');
        
        
        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,V_Psd,'ob','linewidth',1);
        axis([0 100,-0.01,0.5]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('PSD');
        legend('Current-based');
        
%         figure
%         hold on
%         plot(freqbin,V_Psd,'ob','linewidth',1);
%         plot(freqbin,smooth(V_Psd,'moving',15),'-b','linewidth',1);
%         axis([0 100,-0.01,0.1]);
%         set(gca,'fontsize',20);
%         xlabel('Freq.  [Hz]');
%         ylabel('PSD');
        
%         [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vsi);
%         V_Psd = PsdManSmooth;
% 
%         Zi = V_Psd./ZAP_Psd;

%         figure
%         hold on
%         plot(-100,-100,'ob','linewidth',2);
%         plot(-100,-100,'or','linewidth',2);
%         plot(freqbin,Zi,'ob','linewidth',1);
%         %plot(freqbin,Z,'-b','linewidth',1);
%         axis([0 100,-0.01,20]);
%         set(gca,'fontsize',20);
%         xlabel('Freq.  [Hz]');
%         ylabel('Z');
%         
%          figure
%         hold on
%         plot(-100,-100,'ob','linewidth',2);
%         plot(-100,-100,'or','linewidth',2);
%         plot(freqbin,V_Psd,'ob','linewidth',1);
%         %plot(freqbin,Z,'-b','linewidth',1);
%         axis([0 100,-0.01,2]);
%         set(gca,'fontsize',20);
%         xlabel('Freq.  [Hz]');
%         ylabel('Z');

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vscbe);
        V_Psd = PsdManSmooth;
        
        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,V_Psd,'ob','linewidth',1);
        axis([0 100,-0.01,0.5]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('PSD');
        legend('Conductance-based');
        
        
    elseif SYN == 2
        
        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,Isyntrain);
        ZAP_Psd = PsdManSmooth;

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vs);
        Vs_Psd = PsdManSmooth;

        Ze = Vs_Psd./ZAP_Psd;

        figure
        hold on
        plot(freqbin,Ze,'ob','linewidth',1);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,20]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('Z');
        legend('Current-based');
        
        
%         figure
%         hold on
%         plot(-100,-100,'ob','linewidth',2);
%         plot(-100,-100,'or','linewidth',2);
%         plot(freqbin,Vs_Psd,'ob','linewidth',1);
%         axis([0 100,-0.01,0.5]);
%         set(gca,'fontsize',20);
%         xlabel('Freq.  [Hz]');
%         ylabel('PSD');
%         legend('Current-based');
        
        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vscb);
        Vscb_Psd = PsdManSmooth;
        
%         figure
%         hold on
%         plot(-100,-100,'ob','linewidth',2);
%         plot(-100,-100,'or','linewidth',2);
%         plot(freqbin,Vscb_Psd,'ob','linewidth',1);
%         axis([0 100,-0.01,0.5]);
%         set(gca,'fontsize',20);
%         xlabel('Freq.  [Hz]');
%         ylabel('PSD');
%         legend('Conductance-based');
        
        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,Vs_Psd,'ob','linewidth',1);
        plot(freqbin,Vscb_Psd,'or','linewidth',1);
        axis([0 100,-0.01,0.5]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('PSD');
        legend('Current-based','conductance-based');
        
    end
   
    PWSAV = 0;
    if PWSAV == 1
        
        Vs_Psd = zeros(1,length(Vs_Psd));
        Vscb_Psd = zeros(1,length(Vscb_Psd));
        
        for l=1:Ntrials
            
            l
            
            [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vs_vec(l,:));
            Vs_Psd = Vs_Psd+PsdManSmooth;
            
            [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vscb_vec(l,:));
            Vscb_Psd = Vscb_Psd+PsdManSmooth;          
        end
        
        Vs_Psd = Vs_Psd/Ntrials;
        Vscb_Psd = Vscb_Psd/Ntrials;
        
        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,Vs_Psd,'ob','linewidth',1);
        plot(freqbin,Vscb_Psd,'or','linewidth',1);
        axis([0 100,-0.01,0.5]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('PSD');
        legend('Current-based','Conductance-based');
        
    end
   
elseif SYN == 3
    
    [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,Isyntrain);
    ZAP_Psd = PsdManSmooth;

    [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vselin);
    Vlin_Psd = PsdManSmooth;

    Zelin = Vlin_Psd./ZAP_Psd;  

    [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vselin);
    Vpwl_Psd = PsdManSmooth;

    Zepwl = Vpwl_Psd./ZAP_Psd;

    figure
    hold on
%     plot(freqbin,Zepwl,'ob','linewidth',1);
%     plot(freqbin,Zelin,'o','Color',lightcoral,'linewidth',1);
    plot(freqbin,Zepwl,'-b','linewidth',2);
    plot(freqbin,Zelin,'--','Color',lightcoral,'linewidth',2);
    axis([0 100,-0.01,20]);
    set(gca,'fontsize',24);
    xlabel('Freq.  [Hz]');
    ylabel('Z');
    legend('Current-based');

    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(freqbin,Vpwl_Psd,'-b','linewidth',2);
    plot(freqbin,Vlin_Psd,'-','Color',lightcoral,'linewidth',2);
    axis([0 100,-0.01,0.5]);
    set(gca,'fontsize',24);
    xlabel('Freq.  [Hz]');
    ylabel('PSD');
    legend('Current-based');
    
end

% Dependence of the Peak and Troughs values and times as a function of the
% initial conditions

SNG = 0;
if SNG == 1
    
    
    
    
end




