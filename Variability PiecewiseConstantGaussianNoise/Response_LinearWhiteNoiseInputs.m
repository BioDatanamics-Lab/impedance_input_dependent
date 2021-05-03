% 2021-02-17

clearvars;
close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightsalmon = [0.9 0.5 0.4];
mediumacquamarine = [0.4 0.8 0.6];
lightgray = [.7 .7 .7];
darkgray = [.3 .3 .3];
darkgray2 = [.1 .1 .1];

CNT = 0;       
PWC = 1;
        % 1: Linear system
        % 2: Quadratic system (too sensitive)
        % 3: Piecewise linear system

% CNT = 1:  Continuous noise
% PWC = 1, 2, 3:  Piecewise constant noise


% Functions

heaviside=@(t) 0.5*(t == 0)+(t > 0);
pwlv=@(v,gL,gc,alpha) v.*(v<alpha)+(alpha+gc/gL*(v-alpha)).*(v>=alpha);
% gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
gaus = @(x,mu,sig) 1/(sig*sqrt(2*pi))*exp(-(((x-mu).^2)/(2*sig.^2)));


C = 1;
gL = 0.25;
g1 = 0.25;
tau1 = 100;
% 
%     gL = 0.05;
%     g1 = 0.3;
% % 

   g1 = 0;
 
%



% C = 1;
% gL = 0.1;
% g1 = 0.5;
% tau1 = 100;
% 
%   g1 = 0.8;
 %   g1 = 0;

% C = 1;
% gL = 0.1;
% g1 = 5;
% tau1 = 10;

% g1 = 0.8;
% 
% g1 = 0;

% fnat != fres
%  gL = 0.25;
%  g1 = 1.8;
%  tau1 = 100;
 
% fnat != fres
%  gL = 0.3;
%  g1 = 1.5;
%  tau1 = 60;
 
%  gL = 0.3;
%  g1 = 1.3;
%  tau1 = 60;
 
%  gL = 0.3;
%  g1 = 1.3;
%  tau1 = 60;
%  
%  gL = 0.2;
%  g1 = 2;
%  tau1 = 10;

D = 1;

if PWC == 2
    gma = 0;
end

if PWC == 3
    
    alpha = 100;
    gc = 0.1;
    gL = 0.4;
    g1 = 1.3;
    I = 0.5;
    
    tau1 = 60;
    
    % Eigenvalues of the autonomous system 

    a = -gc/C;
    b = -g1/C;
    c = 1/tau1;
    d = -1/tau1;

    [r,mu,fnatc] = Eigenvalues2D(a,b,c,d);

    % Impedance profile

    freq = 1:1:100;
    [Zanl,Phi,fresc,Zmax,QZ] = Impedance2D(a,b,c,d,freq);

    
end

% 0 1 4

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

Tmax = 10000;
dt = 0.1;
t = 0:dt:Tmax;

% Computation of the numerical solution: constant input (or effect of
    % initial conditions

va = zeros(1,length(t));
wa = zeros(1,length(t));

va(1) = 0;
wa(1) = 0;
I = 2;

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

[vamax,jamax] = max(va);
tamax = jamax*dt;
vabar = va(end);
for j=1:length(t);
    if va(j)>=vamax*0.63
        taumrse = t(j);
        break;
    end
end
% if g1==0
%     vamax = vamax*0.8;
%     for j=1:length(t);
%         if va(j)>=vamax
%             tamax = t(j);
%             break
%         end
%     end
% end
% taumrse = vamax/tamax;


figure
hold on
plot(t,va,'-b','linewidth',2);
axis([0 200 0 10]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
ylabel('V');


if PWC == 1
    
    % Time definitions

    Tmax = 1000;
    dt = 0.1;
    t = 0:dt:Tmax;
    
    Tdur = 1;
    Npieces = floor(Tmax/Tdur);
    ton = 0:Tdur:Tmax;
    ton(1) = dt;
    jon = floor(ton/dt);
    
   
    % Random input amplitudes
    
    eta = randn(1,Npieces);
    
    % Uniform input amplitudes 
    
    etamax = 2;
    etamin = -2;
    eta = etamin:4/Npieces:(Npieces-1)*etamax/Npieces;
%     
% %     eta = -2:4/Npieces:(Npieces-1)*2/Npieces;
% %     eta = eta(randperm(length(eta),length(eta)));
% 
%     % Deterministic Gaussian amplitudes
% 
%     etamax = 2;
%     etamin = -2;
%     etaaux = etamin:(etamax-etamin)/(Npieces):(Npieces)*etamax/Npieces;
%     pd = makedist('Normal','mu',0,'sigma',1);
%     eta_cdf = cdf(pd,etaaux); 
%     AmpInt = 2*flip(eta_cdf(1:Npieces/2));
%     AmpInt = AmpInt*(etamax-etamin)/(2*sum(AmpInt));
%     eta(1) = etamin;
%     for j=2:Npieces/2
%         eta(j) = eta(j-1)+AmpInt(j-1);
%     end
%     for j=1:Npieces/2
%         eta(Npieces/2+j) = -eta(Npieces/2-j+1);
%     end
%        
    
    
    Ieta = zeros(1,length(t));    
    for l=1:length(ton)-1        
        Ieta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = eta(l);
    end
    Ieta(end) = Ieta(end-1);
    

    
    % Computation of the solution in response to Ieta
    
    v = zeros(1,length(t));
    w = zeros(1,length(t));

    v(1) = Ieta(1)/(gL+g1);
    w(1) = Ieta(1)/(gL+g1);
    v(1) = v(1)*1.2;

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
    
    veta = v;
    weta = w;
    
    [peakt,vetamax,cntpeak] = PeaksOsc(veta,t,0,Tmax,dt);
    [trought,vetamin,cnttrough] = TroughsOsc(veta,t,0,Tmax,dt);
    cntpeak;
    cnttrough;
    
    venv = zeros(1,Npieces);
    tenv = zeros(1,Npieces);
    for k=2:Npieces
        if eta(k)>=eta(k-1)
            [venv(k),tenv(k)] = max(veta(jon(k):jon(k+1)-1));
            tenv(k) = tenv(k)*dt+ton(k);
        else
            [venv(k),tenv(k)] = min(veta(jon(k):jon(k+1)-1));
            tenv(k) = tenv(k)*dt+ton(k);
        end
    end
   
    
    figure
    hold on
    plot(t,veta,'-b','linewidth',2);
    plot(t,Ieta,'-','Color',lightgray,'linewidth',1);
%     plot(peakt,vetamax,'or','linewidth',2)
%     plot(trought,vetamin,'og','lfinewidth',2)
    plot(tenv,venv,'or','linewidth',2)
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    
    % Ordering of eta according to the increasing amplitude 
    
    etaord = eta;
    L = 0;
    for j=1:Npieces;
        L = L+1;
        for k=L:Npieces
            if etaord(k)<etaord(j)
                aux = etaord(j);
                etaord(j) = etaord(k);
                etaord(k) = aux;
            end
        end 
    end
    
    Ieta = zeros(1,length(t));      
    for l=1:length(ton)-1        
        Ieta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = eta(l);
    end
    Ieta(end) = Ieta(end-1);
    
   
    Ietaord = zeros(1,length(t));      
    for l=1:length(ton)-1        
        Ietaord(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = etaord(l);
    end
    Ietaord(end) = Ietaord(end-1);
    
    % Computation of the solution in response to Ietaord
    
    v = zeros(1,length(t));
    w = zeros(1,length(t));

    v(1) = Ietaord(1)/(gL+g1);
    w(1) = Ietaord(1)/(gL+g1);
    v(1) = v(1)*1.2;

    for j=1:length(t)-1
        k1v = -gL*v(j)-g1*w(j)+D*Ietaord(j);
        k1w = (v(j)-w(j))/tau1;
        av = v(j)+k1v*dt;
        aw = w(j)+k1w*dt;
        k2v = -gL*av-g1*aw+D*Ietaord(j+1);
        k2w = (av-aw)/tau1;
        v(j+1) = v(j)+(k1v+k2v)*dt/2;
        w(j+1) = w(j)+(k1w+k2w)*dt/2;
    end
    
    vord = v;
    word = w;
    
    [peakt,vordmax,cntpeak] = PeaksOsc(vord,t,0,Tmax,dt);
    [trought,vordmin,cnttrough] = TroughsOsc(vord,t,0,Tmax,dt);
    cntpeak;
    cnttrough;
   
    
    if g1 == 0
        vordmax(1) = vord(1);
        vordmax(2:length(ton)) = vord(ton(2:end)/dt);
        peakt = ton(1:end);
    end
    
    vordenv = zeros(1,Npieces);
    tordenv = zeros(1,Npieces);
    for k=2:Npieces
        if etaord(k)>=etaord(k-1)
            [vordenv(k),tordenv(k)] = max(vord(jon(k):jon(k+1)-1));
            tordenv(k) = tordenv(k)*dt+ton(k);
        else
            [vordenv(k),tordenv(k)] = min(vord(jon(k):jon(k+1)-1));
            tordenv(k) = tordenv(k)*dt+ton(k);
        end
    end
   
    figure
    hold on
    plot(t,vord,'-b','linewidth',2);
    plot(t,Ietaord,'-','Color',lightgray,'linewidth',2);
%     plot(peakt,vordmax,'or','linewidth',2)
%     plot(trought,vordmin,'og','linewidth',2)
    plot(tordenv,vordenv,'or','linewidth',2)
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    
    
    SEQ = 0;
            % 1: one input sequence (and all permutations) for each
            %    parameter set
            % 2: Common input sequence (and all permutations)for multiple
            %    parameter sets: changes in gL and g1
            % 3: Common input sequence (and all permutations)for multiple
            %    parameter sets: changes in tau1 and g1
    if SEQ == 1
    
        %   Sequence of random permutations of the order of constant pieces

        Ntrials = 1000;



        vpermsum = zeros(1,length(t));

        vpermenvord_vec = zeros(Ntrials,Npieces);
        tpermenvord_vec = zeros(Ntrials,Npieces);
        for l=1:Ntrials

            l 

            P = randperm(length(etaord),length(etaord));

            etaperm = etaord(P);
            Ietaperm = zeros(1,length(t));    
            for k=1:length(ton)-1        
                Ietaperm(floor(ton(k)/dt):floor((ton(k)+Tdur)/dt)) = etaperm(k);
            end
            Ietaperm(end) = Ietaperm(end-1);

            % Computation of the solution in response to Ietaperm

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end

            if g1 == 0
            end

            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec(l,P) = vpermenv;

            vpermsum = vpermsum+vperm;

    %         jonx = jon;
    %         jonx(1) = 0;
    %         for k=1:length(etaperm)
    %             vpermord(jonx(P(k))+1:jonx(P(k)+1))=vperm(jonx(k)+1:jonx(k+1));
    %         end
            %vpermord(jon(1):jon(2)-1) = vperm(jon(7):jon(8)-1);


        end

        vpermenvordmax = zeros(1,Npieces);
        vpermenvordmin = zeros(1,Npieces);
        vpermenvordmean = zeros(1,Npieces);
        vpermenvordvar = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax(k) = max(vpermenvord_vec(:,k));
            vpermenvordmin(k) = min(vpermenvord_vec(:,k));
            vpermenvordmean(k) = mean(vpermenvord_vec(:,k));
            vpermenvordvar(k) = var(vpermenvord_vec(:,k));
        end

        figure
        hold on
        plot(-100,-100,'-b','linewidth',2);
        %plot(-100,-100,'-r','linewidth',2);
        plot(-100,-100,'-g','linewidth',2);
        for l=1:Ntrials
            plot(vpermenvord_vec(l,:),'o');
        end
        plot(vpermenvordmean,'-b','linewidth',3)
        plot(vordenv(2:end),'-g','linewidth',3) 
        axis([0 Npieces -15 15]);
        set(gca,'fontsize',24);
        xlabel('Piece #');
        ylabel('');
        legend('<P_{\eta}>','P_{\eta,step}')

        figure
        hold on
        plot(t,vperm,'-b','linewidth',2);
        plot(tpermenv,vpermenv,'o','Color',lightcoral,'linewidth',2)
        set(gca,'fontsize',24);
        xlabel('t  [ms]');
        ylabel('V  [mV]');
        legend('V','P_{\eta}');

        figure
        hold on
        plot(vpermenvordvar,'-b','linewidth',2)
        plot(vpermenvordvar/max(va),'-r','linewidth',2)
        set(gca,'fontsize',24);
        xlabel('Piece #');
        ylabel('var(P{\eta})');
        legend('Var(P_{\eta})','VarN(P_{\eta})');

%         figure
%         hold on
%         plot(-100,-100,'-b','linewidth',2);
%         plot(-100,-100,'-r','linewidth',2);
%         plot(-100,-100,'-g','linewidth',2);
%         axis([0 Npieces -15 15]);
%         plot(vpermenvordmax,'-r','linewidth',2)
%         plot(vpermenvordmin,'-r','linewidth',2)
%         plot(vpermenvordmean,'-b','linewidth',2)
%         plot(vordenv(2:end),'-g','linewidth',2) 
%         set(gca,'fontsize',24);
%         xlabel('Piece #');
%         legend('<P_{\eta}>','P_{\eta,ENV}','P_{\eta,step}')

         Ampvpermenvord = (vpermenvordmax-vpermenvordmin)/max(va);

         Data = vpermenvordvar';
         Datanorm = vpermenvordvar'/max(va);
         %save('DataVarRaw_XX_01.dat','Data','-ascii');
         %save('DataVarNorm_XX_01.dat','Datanorm','-ascii');
         
    elseif SEQ == 2
        
        g1vec = [0 2.5 5];
        gLvec = [0.1 0.2];
        
        Ntrials = 1000;



        vpermsum1 = zeros(1,length(t));
        vpermsum2 = zeros(1,length(t));
        vpermsum3 = zeros(1,length(t));
        vpermsum4 = zeros(1,length(t));
        vpermsum5 = zeros(1,length(t));
        vpermsum6 = zeros(1,length(t));

        vpermenvord_vec1 = zeros(Ntrials,Npieces);
        vpermenvord_vec2 = zeros(Ntrials,Npieces);
        vpermenvord_vec3 = zeros(Ntrials,Npieces);
        vpermenvord_vec4 = zeros(Ntrials,Npieces);
        vpermenvord_vec5 = zeros(Ntrials,Npieces);
        vpermenvord_vec6 = zeros(Ntrials,Npieces);
        tpermenvord_vec = zeros(Ntrials,Npieces);
        
        for l=1:Ntrials

            l 

            P = randperm(length(etaord),length(etaord));

            etaperm = etaord(P);
            Ietaperm = zeros(1,length(t));    
            for k=1:length(ton)-1        
                Ietaperm(floor(ton(k)/dt):floor((ton(k)+Tdur)/dt)) = etaperm(k);
            end
            Ietaperm(end) = Ietaperm(end-1);
        
            % Computation of the solution in response to Ietaperm
            
            % 
            gL = gLvec(1);
            g1 = g1vec(1);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec1(l,P) = vpermenv;

            vpermsum1 = vpermsum1+vperm;
            
            % 
            gL = gLvec(1);
            g1 = g1vec(2);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec2(l,P) = vpermenv;

            vpermsum2 = vpermsum1+vperm;

        
            % 
            gL = gLvec(1);
            g1 = g1vec(3);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec3(l,P) = vpermenv;

            vpermsum3 = vpermsum1+vperm;
            
            % 
            gL = gLvec(2);
            g1 = g1vec(1);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec4(l,P) = vpermenv;

            vpermsum4 = vpermsum1+vperm;
            
            % 
            gL = gLvec(2);
            g1 = g1vec(2);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec5(l,P) = vpermenv;

            vpermsum5 = vpermsum1+vperm;
            
            % 
            gL = gLvec(2);
            g1 = g1vec(3);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec6(l,P) = vpermenv;

            vpermsum6 = vpermsum1+vperm;
        
        end
        
        vpermenvordmax1 = zeros(1,Npieces);
        vpermenvordmin1 = zeros(1,Npieces);
        vpermenvordmean1 = zeros(1,Npieces);
        vpermenvordvar1 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax1(k) = max(vpermenvord_vec1(:,k));
            vpermenvordmin1(k) = min(vpermenvord_vec1(:,k));
            vpermenvordmean1(k) = mean(vpermenvord_vec1(:,k));
            vpermenvordvar1(k) = var(vpermenvord_vec1(:,k));
        end
        
        vpermenvordmax2 = zeros(1,Npieces);
        vpermenvordmin2 = zeros(1,Npieces);
        vpermenvordmean2 = zeros(1,Npieces);
        vpermenvordvar2 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax2(k) = max(vpermenvord_vec2(:,k));
            vpermenvordmin2(k) = min(vpermenvord_vec2(:,k));
            vpermenvordmean2(k) = mean(vpermenvord_vec2(:,k));
            vpermenvordvar2(k) = var(vpermenvord_vec2(:,k));
        end
        
        vpermenvordmax3 = zeros(1,Npieces);
        vpermenvordmin3 = zeros(1,Npieces);
        vpermenvordmean3 = zeros(1,Npieces);
        vpermenvordvar3 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax3(k) = max(vpermenvord_vec3(:,k));
            vpermenvordmin3(k) = min(vpermenvord_vec3(:,k));
            vpermenvordmean3(k) = mean(vpermenvord_vec3(:,k));
            vpermenvordvar3(k) = var(vpermenvord_vec3(:,k));
        end
        
        vpermenvordmax4 = zeros(1,Npieces);
        vpermenvordmin4 = zeros(1,Npieces);
        vpermenvordmean4 = zeros(1,Npieces);
        vpermenvordvar4 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax4(k) = max(vpermenvord_vec4(:,k));
            vpermenvordmin4(k) = min(vpermenvord_vec4(:,k));
            vpermenvordmean4(k) = mean(vpermenvord_vec4(:,k));
            vpermenvordvar4(k) = var(vpermenvord_vec4(:,k));
        end
        
        vpermenvordmax5 = zeros(1,Npieces);
        vpermenvordmin5 = zeros(1,Npieces);
        vpermenvordmean5 = zeros(1,Npieces);
        vpermenvordvar5 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax5(k) = max(vpermenvord_vec5(:,k));
            vpermenvordmin5(k) = min(vpermenvord_vec5(:,k));
            vpermenvordmean5(k) = mean(vpermenvord_vec5(:,k));
            vpermenvordvar5(k) = var(vpermenvord_vec5(:,k));
        end
        
        vpermenvordmax6 = zeros(1,Npieces);
        vpermenvordmin6 = zeros(1,Npieces);
        vpermenvordmean6 = zeros(1,Npieces);
        vpermenvordvar6 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax6(k) = max(vpermenvord_vec6(:,k));
            vpermenvordmin6(k) = min(vpermenvord_vec6(:,k));
            vpermenvordmean6(k) = mean(vpermenvord_vec6(:,k));
            vpermenvordvar6(k) = var(vpermenvord_vec6(:,k));
        end
        
        figure
        hold on
        plot(vpermenvordvar1,'-b','linewidth',2)
        plot(vpermenvordvar2,'-r','linewidth',2)
        plot(vpermenvordvar3,'-g','linewidth',2)
        set(gca,'fontsize',24);
        xlabel('Piece #');
        ylabel('var(P{\eta})');
        legend('g_1=0','g_1=2.5','g_1=5');
        title('g_L = 0.1');
        
        figure
        hold on
        plot(vpermenvordvar4,'-b','linewidth',2)
        plot(vpermenvordvar5,'-r','linewidth',2)
        plot(vpermenvordvar6,'-g','linewidth',2)
        set(gca,'fontsize',24);
        xlabel('Piece #');
        ylabel('var(P{\eta})');
        legend('g_1=0','g_1=2.5','g_1=5');
        title('g_L = 0.2');
        
        figure
        hold on
        plot(vpermenvordvar1,'ob','linewidth',2)
        plot(vpermenvordvar2,'or','linewidth',2)
        plot(vpermenvordvar3,'og','linewidth',2)
        plot(vpermenvordvar4,'ob','linewidth',2)
        plot(vpermenvordvar5,'or','linewidth',2)
        plot(vpermenvordvar6,'og','linewidth',2)
        set(gca,'fontsize',24);
        xlabel('Piece #');
        ylabel('var(P{\eta})');
        legend('g1=0','g1=2.5','g1=5');
       


        Data1 = vpermenvordvar1';
        Data2 = vpermenvordvar2';
        Data3 = vpermenvordvar3';
        Data4 = vpermenvordvar4';
        Data5 = vpermenvordvar5';
        Data6 = vpermenvordvar6';
           
    elseif SEQ == 3
        
        gL = 0.4;
        g1vec = [1 2];
        tau1vec = [100 50 25];
        
        Ntrials = 1000;



        vpermsum1 = zeros(1,length(t));
        vpermsum2 = zeros(1,length(t));
        vpermsum3 = zeros(1,length(t));
        vpermsum4 = zeros(1,length(t));
        vpermsum5 = zeros(1,length(t));
        vpermsum6 = zeros(1,length(t));

        vpermenvord_vec1 = zeros(Ntrials,Npieces);
        vpermenvord_vec2 = zeros(Ntrials,Npieces);
        vpermenvord_vec3 = zeros(Ntrials,Npieces);
        vpermenvord_vec4 = zeros(Ntrials,Npieces);
        vpermenvord_vec5 = zeros(Ntrials,Npieces);
        vpermenvord_vec6 = zeros(Ntrials,Npieces);
        tpermenvord_vec = zeros(Ntrials,Npieces);
        
        for l=1:Ntrials

            l 

            P = randperm(length(etaord),length(etaord));

            etaperm = etaord(P);
            Ietaperm = zeros(1,length(t));    
            for k=1:length(ton)-1        
                Ietaperm(floor(ton(k)/dt):floor((ton(k)+Tdur)/dt)) = etaperm(k);
            end
            Ietaperm(end) = Ietaperm(end-1);
        
            % Computation of the solution in response to Ietaperm
            
            % 
            
            g1 = g1vec(1);
            tau1 = tau1vec(1);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec1(l,P) = vpermenv;

            vpermsum1 = vpermsum1+vperm;
            
            %            
            g1 = g1vec(1);
            tau1 = tau1vec(2);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec2(l,P) = vpermenv;

            vpermsum2 = vpermsum1+vperm;

        
            % 
            g1 = g1vec(1);
            tau1 = tau1vec(3);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec3(l,P) = vpermenv;

            vpermsum3 = vpermsum1+vperm;
            
            % 
            g1 = g1vec(2);
            tau1 = tau1vec(1);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec4(l,P) = vpermenv;

            vpermsum4 = vpermsum1+vperm;
            
            % 
            g1 = g1vec(2);
            tau1 = tau1vec(2);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec5(l,P) = vpermenv;

            vpermsum5 = vpermsum1+vperm;
            
            % 
            g1 = g1vec(2);
            tau1 = tau1vec(3);

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end
            
            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end

            vpermenvord_vec6(l,P) = vpermenv;

            vpermsum6 = vpermsum1+vperm;
        
        end
        
        vpermenvordmax1 = zeros(1,Npieces);
        vpermenvordmin1 = zeros(1,Npieces);
        vpermenvordmean1 = zeros(1,Npieces);
        vpermenvordvar1 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax1(k) = max(vpermenvord_vec1(:,k));
            vpermenvordmin1(k) = min(vpermenvord_vec1(:,k));
            vpermenvordmean1(k) = mean(vpermenvord_vec1(:,k));
            vpermenvordvar1(k) = var(vpermenvord_vec1(:,k));
        end
        
        vpermenvordmax2 = zeros(1,Npieces);
        vpermenvordmin2 = zeros(1,Npieces);
        vpermenvordmean2 = zeros(1,Npieces);
        vpermenvordvar2 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax2(k) = max(vpermenvord_vec2(:,k));
            vpermenvordmin2(k) = min(vpermenvord_vec2(:,k));
            vpermenvordmean2(k) = mean(vpermenvord_vec2(:,k));
            vpermenvordvar2(k) = var(vpermenvord_vec2(:,k));
        end
        
        vpermenvordmax3 = zeros(1,Npieces);
        vpermenvordmin3 = zeros(1,Npieces);
        vpermenvordmean3 = zeros(1,Npieces);
        vpermenvordvar3 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax3(k) = max(vpermenvord_vec3(:,k));
            vpermenvordmin3(k) = min(vpermenvord_vec3(:,k));
            vpermenvordmean3(k) = mean(vpermenvord_vec3(:,k));
            vpermenvordvar3(k) = var(vpermenvord_vec3(:,k));
        end
        
        vpermenvordmax4 = zeros(1,Npieces);
        vpermenvordmin4 = zeros(1,Npieces);
        vpermenvordmean4 = zeros(1,Npieces);
        vpermenvordvar4 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax4(k) = max(vpermenvord_vec4(:,k));
            vpermenvordmin4(k) = min(vpermenvord_vec4(:,k));
            vpermenvordmean4(k) = mean(vpermenvord_vec4(:,k));
            vpermenvordvar4(k) = var(vpermenvord_vec4(:,k));
        end
        
        vpermenvordmax5 = zeros(1,Npieces);
        vpermenvordmin5 = zeros(1,Npieces);
        vpermenvordmean5 = zeros(1,Npieces);
        vpermenvordvar5 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax5(k) = max(vpermenvord_vec5(:,k));
            vpermenvordmin5(k) = min(vpermenvord_vec5(:,k));
            vpermenvordmean5(k) = mean(vpermenvord_vec5(:,k));
            vpermenvordvar5(k) = var(vpermenvord_vec5(:,k));
        end
        
        vpermenvordmax6 = zeros(1,Npieces);
        vpermenvordmin6 = zeros(1,Npieces);
        vpermenvordmean6 = zeros(1,Npieces);
        vpermenvordvar6 = zeros(1,Npieces);
        for k=1:Npieces
            vpermenvordmax6(k) = max(vpermenvord_vec6(:,k));
            vpermenvordmin6(k) = min(vpermenvord_vec6(:,k));
            vpermenvordmean6(k) = mean(vpermenvord_vec6(:,k));
            vpermenvordvar6(k) = var(vpermenvord_vec6(:,k));
        end
        
        figure
        hold on
        plot(vpermenvordvar1,'-b','linewidth',2)
        plot(vpermenvordvar2,'-r','linewidth',2)
        plot(vpermenvordvar3,'-g','linewidth',2)
        set(gca,'fontsize',24);
        xlabel('Piece #');
        ylabel('var(P{\eta})');
        legend('\tau_1=100','\tau_1=50','\tau_1=25');
        title('g_1 = 1');
        
        figure
        hold on
        plot(vpermenvordvar4,'-b','linewidth',2)
        plot(vpermenvordvar5,'-r','linewidth',2)
        plot(vpermenvordvar6,'-g','linewidth',2)
        set(gca,'fontsize',24);
        xlabel('Piece #');
        ylabel('var(P{\eta})');
        legend('\tau_1=100','\tau_1=50','\tau_1=25');
        title('g_1 = 2');
        
        figure
        hold on
        plot(vpermenvordvar1,'ob','linewidth',2)
        plot(vpermenvordvar2,'or','linewidth',2)
        plot(vpermenvordvar3,'og','linewidth',2)
        plot(vpermenvordvar4,'ob','linewidth',2)
        plot(vpermenvordvar5,'or','linewidth',2)
        plot(vpermenvordvar6,'og','linewidth',2)
        set(gca,'fontsize',24);
        xlabel('Piece #');
        ylabel('var(P{\eta})');
        legend('\tau_1=100','\tau_1=50','\tau_1=25');

       


        Data1 = vpermenvordvar1';
        Data2 = vpermenvordvar2';
        Data3 = vpermenvordvar3';
        Data4 = vpermenvordvar4';
        Data5 = vpermenvordvar5';
        Data6 = vpermenvordvar6';
        
        
        
    end
    
    VRNC = 0;
            % 1: Single values of g1 & gL
            % 2: Effects of changes in g1
            % 3: Effects of changes in gL
    if VRNC == 1
        
        Ntrials = 100;

        for l=1:Ntrials

            l 

            P = randperm(length(etaord),length(etaord));

            etaperm = etaord(P);
            Ietaperm = zeros(1,length(t));    
            for k=1:length(ton)-1        
                Ietaperm(floor(ton(k)/dt):floor((ton(k)+Tdur)/dt)) = etaperm(k);
            end
            Ietaperm(end) = Ietaperm(end-1);

            % Computation of the solution in response to Ietaperm

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = Ietaperm(1)/(gL+g1);
            w(1) = Ietaperm(1)/(gL+g1);
            v(1) = v(1)*1.2;

            for j=1:length(t)-1
                k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            meanv(l) = mean(v);
            varv(l) = var(v);   
            
        end
        
        [mean(meanv) mean(varv)]
            
    elseif VRNC == 2
            
        g1vec = 0:0.1:2; 
        Ntrials = 100;
        
        Ietaperm_vec = zeros(Ntrials,length(etaord));
        
        for l=1:Ntrials

            P = randperm(length(etaord),length(etaord));

            etaperm = etaord(P);
            Ietaperm = zeros(1,length(t));    
            for k=1:length(ton)-1        
                Ietaperm(floor(ton(k)/dt):floor((ton(k)+Tdur)/dt)) = etaperm(k);
            end
            Ietaperm(end) = Ietaperm(end-1);
            
            Ietaperm_vec(l,1:length(Ietaperm)) = Ietaperm;
        end
        
        v = zeros(1,length(t));
        w = zeros(1,length(t));
        mean_vec = zeros(1,length(g1vec));
        var_vec = zeros(1,length(g1vec));
        meanv = zeros(1,Ntrials);
        varv = zeros(1,Ntrials);
        
        for i=1:length(g1vec)
        
            g1 = g1vec(i);
            
            for l=1:Ntrials     

                v(1) = Ietaperm(1)/(gL+g1);
                w(1) = Ietaperm(1)/(gL+g1);
                v(1) = v(1)*1.2;

                for j=1:length(t)-1
                    k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                    k1w = (v(j)-w(j))/tau1;
                    av = v(j)+k1v*dt;
                    aw = w(j)+k1w*dt;
                    k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                    k2w = (av-aw)/tau1;
                    v(j+1) = v(j)+(k1v+k2v)*dt/2;
                    w(j+1) = w(j)+(k1w+k2w)*dt/2;
                end

                meanv(l) = mean(v);
                varv(l) = var(v);   
            end
            
            mean_vec(i) = mean(meanv);
            var_vec(i) = mean(varv);
            
            [g1vec(i) mean_vec(i) var_vec(i)];
            
        end
        
        figure
        hold on
        plot(g1vec,var_vec,'-o')
        set(gca,'fontsize',24);
        xlabel('g_1');
        ylabel('Var');
             
        
    elseif VRNC == 3
            
        gLvec = 0.1:0.1:2; 
        Ntrials = 100;
        
        Ietaperm_vec = zeros(Ntrials,length(etaord));
        
        for l=1:Ntrials

            P = randperm(length(etaord),length(etaord));

            etaperm = etaord(P);
            Ietaperm = zeros(1,length(t));    
            for k=1:length(ton)-1        
                Ietaperm(floor(ton(k)/dt):floor((ton(k)+Tdur)/dt)) = etaperm(k);
            end
            Ietaperm(end) = Ietaperm(end-1);
            
            Ietaperm_vec(l,1:length(Ietaperm)) = Ietaperm;
        end
        
        v = zeros(1,length(t));
        w = zeros(1,length(t));
        mean_vec = zeros(1,length(gLvec));
        var_vec = zeros(1,length(gLvec));
        meanv = zeros(1,Ntrials);
        varv = zeros(1,Ntrials);
        
        for i=1:length(gLvec)
        
            gL = gLvec(i);
            
            for l=1:Ntrials     

                v(1) = Ietaperm(1)/(gL+g1);
                w(1) = Ietaperm(1)/(gL+g1);
                v(1) = v(1)*1.2;

                for j=1:length(t)-1
                    k1v = -gL*v(j)-g1*w(j)+D*Ietaperm(j);
                    k1w = (v(j)-w(j))/tau1;
                    av = v(j)+k1v*dt;
                    aw = w(j)+k1w*dt;
                    k2v = -gL*av-g1*aw+D*Ietaperm(j+1);
                    k2w = (av-aw)/tau1;
                    v(j+1) = v(j)+(k1v+k2v)*dt/2;
                    w(j+1) = w(j)+(k1w+k2w)*dt/2;
                end

                meanv(l) = mean(v);
                varv(l) = var(v);   
            end
            
            mean_vec(i) = mean(meanv);
            var_vec(i) = mean(varv);
            
            [gLvec(i) mean_vec(i) var_vec(i)]
            
        end
        
        figure
        hold on
        plot(gLvec,var_vec,'-o')
        set(gca,'fontsize',24);
        xlabel('g_1');
        ylabel('Var');
        
    end
 
    PWS = 2;
            % 1: PSD and Z for a PWC signal with random amplitudes
            % 2: PSD and Z for a PWC signal with non-equal, but non-random
            %    amplitudes
            %
            % PWC: Piecewise constant
            % PSD: Power spectrum density 
            % Z:   Impedance
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


        % Random input amplitudes

        eta = randn(1,Npieces);
        
        Ieta = zeros(1,length(t));    
        for l=1:length(ton)-1        
            Ieta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = eta(l);
        end
        Ieta(end) = Ieta(end-1);


        % Computation of the solution in response to Ieta

        v = zeros(1,length(t));
        w = zeros(1,length(t));

        v(1) = Ieta(1)/(gL+g1);
        w(1) = Ieta(1)/(gL+g1);
        v(1) = v(1)*1.2;

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

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,D*Ieta);
        ZAP_Psd = PsdManSmooth;

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,v);
        V_Psd = PsdManSmooth;

        Z = V_Psd./ZAP_Psd;

        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,Z,'ob','linewidth',1);
        plot(freqbin,V_Psd*max(Z)/max(V_Psd),'or','linewidth',1);
        plot(fres,0,'ok','linewidth',3);
        plot(fnat,0,'ok','linewidth',3);
        plot([fres fres],[0 20],'--k','linewidth',2);
        plot([fnat fnat],[0 20],'--k','linewidth',2);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,6]);
        set(gca,'fontsize',20);
        xlabel('Freq.  [Hz]');
        ylabel('Z');
        
        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,V_Psd,'ob','linewidth',1);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,2]);
        set(gca,'fontsize',20);
        xlabel('Freq.  [Hz]');
        ylabel('PSD');



        
    elseif PWS == 2 
        
        

        
        % Time definitions

        Tmaxx = 1000000;
        dt = 0.1;
        t = 0:dt:Tmaxx;

        Tdur = 1;
        Npieces = floor(Tmaxx/Tdur);
        ton = 0:Tdur:Tmaxx;
        ton(1) = dt;
        jon = floor(ton/dt);
        
         % Increasing, non-random input amplitudes
        
        etamax = 2;
        etamin = -2;
        etaord = etamin:4/Npieces:(Npieces-1)*etamax/Npieces;
        
%         % Deterministic Gaussian-like amplitudes
% 
%         etamax = 2;
%         etamin = -2;
%         etaaux = etamin:(etamax-etamin)/(Npieces):(Npieces)*etamax/Npieces;
%         pd = makedist('Normal','mu',0,'sigma',1.5);
%         eta_cdf = cdf(pd,etaaux); 
%         AmpInt = 2*flip(eta_cdf(1:Npieces/2));
%         AmpInt = AmpInt*(etamax-etamin)/(2*sum(AmpInt));
%         eta(1) = etamin;
%         for j=2:Npieces/2
%             eta(j) = eta(j-1)+AmpInt(j-1);
%         end
%         for j=1:Npieces/2
%             eta(Npieces/2+j) = -eta(Npieces/2-j+1);
%         end
%         etaord = eta;
        
        Ietaord = zeros(1,length(t));    
        for l=1:length(ton)-1        
            Ietaord(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = etaord(l);
        end
        Ietaord(end) = Ietaord(end-1);
        
        
        v = zeros(1,length(t));
        w = zeros(1,length(t));

        v(1) = Ietaord(1)/(gL+g1);
        w(1) = Ietaord(1)/(gL+g1);
        

        for j=1:length(t)-1
            k1v = -gL*v(j)-g1*w(j)+D*Ietaord(j);
            k1w = (v(j)-w(j))/tau1;
            av = v(j)+k1v*dt;
            aw = w(j)+k1w*dt;
            k2v = -gL*av-g1*aw+D*Ietaord(j+1);
            k2w = (av-aw)/tau1;
            v(j+1) = v(j)+(k1v+k2v)*dt/2;
            w(j+1) = w(j)+(k1w+k2w)*dt/2;
        end
        
        vord = v;
        word = w;
        
        venv = zeros(1,Npieces);
        tenv = zeros(1,Npieces);
        for k=2:Npieces
            if etaord(k)>=etaord(k-1)
                [venv(k),tenv(k)] = max(vord(jon(k):jon(k+1)-1));
                tenv(k) = tenv(k)*dt+ton(k);
            else
                [venv(k),tenv(k)] = min(vord(jon(k):jon(k+1)-1));
                tenv(k) = tenv(k)*dt+ton(k);
            end
        end
        
        vordenv = venv;
        tordenv = tenv;

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,D*Ietaord);
        ZAPord_Psd = PsdManSmooth;
        PSD;
        Freq;
        freqbin;

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vord);
        Vord_Psd = PsdManSmooth;
        PSD;
        Freq;
        freqbin;

        Zord = Vord_Psd./ZAPord_Psd;

        
       
        % Input amplitudes in random order (same set as i the previous
        % step)
       
        
        etarand = etaord(randperm(length(etaord),length(etaord)));
        
        Ietarand = zeros(1,length(t));    
        for l=1:length(ton)-1        
            Ietarand(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = etarand(l);
        end
        Ietarand(end) = Ietarand(end-1);
        
        v = zeros(1,length(t));
        w = zeros(1,length(t));

        v(1) = Ietarand(1)/(gL+g1);
        w(1) = Ieta(1)/(gL+g1);

        for j=1:length(t)-1
            k1v = -gL*v(j)-g1*w(j)+D*Ietarand(j);
            k1w = (v(j)-w(j))/tau1;
            av = v(j)+k1v*dt;
            aw = w(j)+k1w*dt;
            k2v = -gL*av-g1*aw+D*Ietarand(j+1);
            k2w = (av-aw)/tau1;
            v(j+1) = v(j)+(k1v+k2v)*dt/2;
            w(j+1) = w(j)+(k1w+k2w)*dt/2;
        end
        
        vrand = v;
        wrand = w;
        
        venv = zeros(1,Npieces);
        tenv = zeros(1,Npieces);
        for k=2:Npieces
            if etarand(k)>=etarand(k-1)
                [venv(k),tenv(k)] = max(vrand(jon(k):jon(k+1)-1));
                tenv(k) = tenv(k)*dt+ton(k);
            else
                [venv(k),tenv(k)] = min(vrand(jon(k):jon(k+1)-1));
                tenv(k) = tenv(k)*dt+ton(k);
            end
        end
        
        vrandenv = venv;
        trandenv = tenv;

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,Ietarand);
        ZAPrand_Psd = PsdManSmooth;
        PSD;
        Freq;
        freqbin;

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vrand);
        Vrand_Psd = PsdManSmooth;
        PSD;
        Freq;
        freqbin;

        Zrand = Vrand_Psd./ZAPrand_Psd;
        
        figure
        hold on
        plot(t,vord,'-b','linewidth',2);
        plot(t,Ietaord,'-','Color',lightgray,'linewidth',1);
        plot(tordenv,vordenv,'or','linewidth',2)
        set(gca,'fontsize',24);
        xlabel('t  [ms]');
        ylabel('V  [mV]');
        
        figure
        hold on
        plot(t,vrand,'-b','linewidth',2);
        plot(t,Ietarand,'-','Color',lightgray,'linewidth',1);
        plot(trandenv,vrandenv,'or','linewidth',2)
        set(gca,'fontsize',24);
        xlabel('t  [ms]');
        ylabel('V  [mV]');
        
        
        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,Zord,'ob','linewidth',1);
        plot(freqbin,Vord_Psd*max(Zord)/max(Vord_Psd),'or','linewidth',1);
        plot(fres,0,'ok','linewidth',3);
        plot(fnat,0,'ok','linewidth',3);
        plot([fres fres],[0 20],'--k','linewidth',2);
        plot([fnat fnat],[0 20],'--k','linewidth',2);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,10]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('Z');
        legend('Z','PSD');
        
        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,Zrand,'ob','linewidth',1);
        plot(freqbin,Vrand_Psd*max(Zrand)/max(Vrand_Psd),'or','linewidth',1);
        plot(fres,0,'ok','linewidth',3);
        plot(fnat,0,'ok','linewidth',3);
        plot([fres fres],[0 20],'--k','linewidth',2);
        plot([fnat fnat],[0 20],'--k','linewidth',2);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,10]);
        set(gca,'fontsize',24);
        xlabel('Freq.  [Hz]');
        ylabel('Z');
        legend('Z','PSD');
        
        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,Zrand,'ob','linewidth',1);
%         plot(fres,0,'ok','linewidth',3);
%         plot(fnat,0,'ok','linewidth',3);
%         plot([fres fres],[0 20],'--k','linewidth',2);
%         plot([fnat fnat],[0 20],'--k','linewidth',2);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,5]);
        set(gca,'fontsize',20);
        xlabel('Freq.  [Hz]');
        ylabel('Z');
        
        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,Vrand_Psd,'ob','linewidth',1);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,0.6]);
        set(gca,'fontsize',20);
        xlabel('Freq.  [Hz]');
        ylabel('PSD');
        
        
%         figure
%         hold on
%         plot(-100,-100,'ob','linewidth',2);
%         plot(-100,-100,'or','linewidth',2);
%         plot(freqbin,Zrand,'ob','linewidth',1);
%         plot(freqbin,Vrand_Psd*max(Zrand)/max(Vrand_Psd),'or','linewidth',1);
%         plot(freqbin,Zord,'o','Color',lightblueish,'linewidth',1);
%         plot(freqbin,Vord_Psd*max(Zord)/max(Vord_Psd),'o','Color',lightcoral,'linewidth',1);
%         plot(fres,0,'ok','linewidth',3);
%         plot(fnat,0,'ok','linewidth',3);
%         plot([fres fres],[0 20],'--k','linewidth',2);
%         plot([fnat fnat],[0 20],'--k','linewidth',2);
%         %plot(freqbin,Z,'-b','linewidth',1);
%         axis([0 100,-0.01,10]);
%         set(gca,'fontsize',24);
%         xlabel('Freq.  [Hz]');
%         ylabel('Z');
%         legend('Z         random order','PSD    random order','Z         increasing order','PSD    increasing order');
%       
%         [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,vord(10000:end));
%         Vordcut_Psd = PsdManSmooth;
%         PSD;
%         Freq;
%         freqbin;
%         plot(freqbin,Vordcut_Psd*max(Zord)/max(Vordcut_Psd),'o','Color',mediumacquamarine,'linewidth',1);
%         legend('Z         random','PSD    random','Z         ordered','PSD    ordered','PSD    ordered (cut)');
%         
%         
        
    end
  
   
    
    
%     figure(2)
%     plot(t,vord,'--b','linewidth',2);
%     plot(t,Ietaord,'--b','Color',lightgray,'linewidth',2);

%     % Constant input amplitudes
%     
%     xeta = min(eta):(max(eta)-min(eta))/(Npieces-1):max(eta);
%     
%     
%     Ixeta = zeros(1,length(t));      
%     for l=1:length(ton)-1        
%         Ixeta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = xeta(l);
%     end
%     Ixeta(end) = Ixeta(end-1);
%     
%     % Computation of the solution in response to Ixeta
%     
%     v = zeros(1,length(t));
%     w = zeros(1,length(t));
% 
%     v(1) = Ixeta(1)/(gL+g1);
%     w(1) = Ixeta(1)/(gL+g1);
%     v(1) = v(1)*1.2;
% 
%     for j=1:length(t)-1
%         k1v = -gL*v(j)-g1*w(j)+D*Ixeta(j);
%         k1w = (v(j)-w(j))/tau1;
%         av = v(j)+k1v*dt;
%         aw = w(j)+k1w*dt;
%         k2v = -gL*av-g1*aw+D*Ixeta(j+1);
%         k2w = (av-aw)/tau1;
%         v(j+1) = v(j)+(k1v+k2v)*dt/2;
%         w(j+1) = w(j)+(k1w+k2w)*dt/2;
%     end
%     
%     vstep = v;
%     wstep = w;
%     
%     [peakt,vstepmax,cntpeak] = PeaksOsc(vstep,t,0,Tmax,dt);
%     [trought,vstepmin,cnttrough] = TroughsOsc(vstep,t,0,Tmax,dt);
%     cntpeak;
%     cnttrough;
%     
%     if g1 == 0
%         vstepmax(1) = vstep(1);
%         vstepmax(2:length(ton)) = vstep(ton(2:end)/dt);
%         peakt = ton(1:end);
%     end
%     
%     figure
%     hold on
%     plot(t,vstep,'-b','linewidth',2);
%     plot(t,Ixeta,'--','Color',lightgray,'linewidth',2);
%     plot(peakt,vstepmax,'or','linewidth',2)
%     plot(trought,vstepmin,'og','linewidth',2)
%     set(gca,'fontsize',24);
%     xlabel('t  [ms]');
%     ylabel('V  [mV]');
    

elseif PWC == 2
    
    gma = 1.4;
    D = 0.1;
    I = 0;
    
    g1 = 0.2;
    
    % Fixed-points
    
    vfpaux = roots([gL*gma -(gL+g1) I]);
    vfp = min(vfpaux);
    wfp = vfp;
    
    % Computation of the numerical solution: constant input (or effect of
    % initial conditions

    va = zeros(1,length(t));
    wa = zeros(1,length(t));

    va(1) = vfp;
    wa(1) = wfp;
    

    for j=1:length(t)-1
        k1v = -gL*va(j)+gL*gma*va(j)^2-g1*wa(j)+I;
        k1w = (va(j)-wa(j))/tau1;
        av = va(j)+k1v*dt;
        aw = wa(j)+k1w*dt;
        k2v = -gL*av+gL*gma*av^2-g1*aw+I;
        k2w = (av-aw)/tau1;
        va(j+1) = va(j)+(k1v+k2v)*dt/2;
        wa(j+1) = wa(j)+(k1w+k2w)*dt/2;
    end
    
    [vamax,jamax] = max(va);
    tamax = jamax*dt;
    vabar = va(end);
    for j=1:length(t);
        if va(j)>=vamax*0.63
            taumrse = t(j);
            break;
        end
    end
    
    figure
    hold on
    plot(t,va,'-b','linewidth',2);
    axis([0 200 -10 10]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    
    vv = -20:0.01:20;
    figure
    hold on
    plot(vv,(-gL*vv+gL*gma*vv.^2+I)/g1,'-r','linewidth',2);
    plot(vv,vv,'-g','linewidth',2);
    plot(va,wa,'-b','linewidth',2);
    plot([-20 20],[0 0],'--k')
    plot([0 0],[-20 20],'--k')
    axis([-5 5 -1 4]);
    set(gca,'fontsize',24);
    xlabel('V  [mV]');
    ylabel('w');
    
    
    % Time definitions

    Tmax = 1000;
    dt = 0.1;
    t = 0:dt:Tmax;
    
    Tdur = 1;
    Npieces = floor(Tmax/Tdur);
    ton = 0:Tdur:Tmax;
    ton(1) = dt;
    jon = floor(ton/dt);
    
   
    % Random input amplitudes
    
    eta = randn(1,Npieces);
    
%     eta = -2:4/Npieces:(Npieces-1)*2/Npieces;
%     eta = eta(randperm(length(eta),length(eta)));
    
    Ieta = zeros(1,length(t));    
    for l=1:length(ton)-1        
        Ieta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = eta(l);
    end
    Ieta(end) = Ieta(end-1);
    
    % Computation of the solution in response to Ieta
    
    v = zeros(1,length(t));
    w = zeros(1,length(t));

    v(1) = vfp;
    w(1) = wfp;

    for j=1:length(t)-1
        k1v = -gL*v(j)+gL*gma*v(j)^2-g1*w(j)+D*Ieta(j);
        k1w = (v(j)-w(j))/tau1;
        av = v(j)+k1v*dt;
        aw = w(j)+k1w*dt;
        k2v = -gL*av+gL*gma*av^2-g1*aw+D*Ieta(j+1);
        k2w = (av-aw)/tau1;
        v(j+1) = v(j)+(k1v+k2v)*dt/2;
        w(j+1) = w(j)+(k1w+k2w)*dt/2;
    end
    
    veta = v;
    weta = w;
    
    venv = zeros(1,Npieces);
    tenv = zeros(1,Npieces);
    for k=2:Npieces
        if eta(k)>=eta(k-1)
            [venv(k),tenv(k)] = max(veta(jon(k):jon(k+1)-1));
            tenv(k) = tenv(k)*dt+ton(k);
        else
            [venv(k),tenv(k)] = min(veta(jon(k):jon(k+1)-1));
            tenv(k) = tenv(k)*dt+ton(k);
        end
    end
   
    figure
    hold on
    plot(t,veta,'-b','linewidth',2);
    plot(t,Ieta,'-','Color',lightgray,'linewidth',1);
    plot(tenv,venv,'or','linewidth',1)
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    
    % Ordering of eta according to the increasing amplitude 
    
    etaord = eta;
    L = 0;
    for j=1:Npieces;
        L = L+1;
        for k=L:Npieces
            if etaord(k)<etaord(j)
                aux = etaord(j);
                etaord(j) = etaord(k);
                etaord(k) = aux;
            end
        end 
    end
    
    Ieta = zeros(1,length(t));      
    for l=1:length(ton)-1        
        Ieta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = eta(l);
    end
    Ieta(end) = Ieta(end-1);
    
   
    Ietaord = zeros(1,length(t));      
    for l=1:length(ton)-1        
        Ietaord(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = etaord(l);
    end
    Ietaord(end) = Ietaord(end-1);
    
    % Computation of the solution in response to Ietaord
    
    v = zeros(1,length(t));
    w = zeros(1,length(t));

    v(1) = vfp;
    w(1) = wfp;

    for j=1:length(t)-1
        k1v = -gL*v(j)+gL*gma*v(j)^2-g1*w(j)+D*Ietaord(j);
        k1w = (v(j)-w(j))/tau1;
        av = v(j)+k1v*dt;
        aw = w(j)+k1w*dt;
        k2v = -gL*av+gL*gma*av^2-g1*aw+D*Ietaord(j+1);
        k2w = (av-aw)/tau1;
        v(j+1) = v(j)+(k1v+k2v)*dt/2;
        w(j+1) = w(j)+(k1w+k2w)*dt/2;
    end
    
    vord = v;
    word = w;
    
    vordenv = zeros(1,Npieces);
    tordenv = zeros(1,Npieces);
    for k=2:Npieces
        if etaord(k)>=etaord(k-1)
            [vordenv(k),tordenv(k)] = max(vord(jon(k):jon(k+1)-1));
            tordenv(k) = tordenv(k)*dt+ton(k);
        else
            [vordenv(k),tordenv(k)] = min(vord(jon(k):jon(k+1)-1));
            tordenv(k) = tordenv(k)*dt+ton(k);
        end
    end
   
    figure
    hold on
    plot(t,vord,'-b','linewidth',2);
    plot(t,Ietaord,'-','Color',lightgray,'linewidth',1);
    plot(tordenv,vordenv,'or','linewidth',1);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    
elseif PWC == 3
    
    
 
    
    pwlv=@(v,gL,gc,alpha) v.*(v<alpha)+(alpha+gc/gL*(v-alpha)).*(v>=alpha);
    
    
    
    % Computation of the numerical solution: constant input (or effect of
    % initial conditions

    va = zeros(1,length(t));
    wa = zeros(1,length(t));

    va(1) = 0.01;
    wa(1) = 0.01;
    

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
    
    
    [vamax,jamax] = max(va);
    tamax = jamax*dt;
    vabar = va(end);
    for j=1:length(t);
        if va(j)>=vamax*0.63
            taumrse = t(j);
            break;
        end
    end
    
    figure
    hold on
    plot(t,va,'-b','linewidth',2);
    axis([0 200 -10 10]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    
    vv = -4:0.01:4;
    figure
    hold on
    plot(vv,-gL*pwlv(vv,gL,gc,alpha)/g1,'-r','linewidth',2);
    plot(vv,vv,'-g','linewidth',2);
    plot(va,wa,'-b','linewidth',2);
    plot(vv,(-gL*pwlv(vv,gL,gc,alpha)+I)/g1,'--r','linewidth',2);
    plot([-20 20],[0 0],'--k')
    plot([0 0],[-20 20],'--k')
    axis([-3 3 -3 3]);
    set(gca,'fontsize',24);
    xlabel('V  [mV]');
    ylabel('w');
    legend('v-nullclne','w-nullcline','trajectory','Location','SouthWest');
    
     % Time definitions

    Tmax = 1000;
    dt = 0.1;
    t = 0:dt:Tmax;
    
    Tdur = 1;
    Npieces = floor(Tmax/Tdur);
    ton = 0:Tdur:Tmax;
    ton(1) = dt;
    jon = floor(ton/dt);
    
   
    % Random input amplitudes
    
    eta = randn(1,Npieces);
    
%     eta = -2:4/Npieces:(Npieces-1)*2/Npieces;
%     eta = eta(randperm(length(eta),length(eta)));
    
    Ieta = zeros(1,length(t));    
    for l=1:length(ton)-1        
        Ieta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = eta(l);
    end
    Ieta(end) = Ieta(end-1);
    
    % Computation of the solution in response to Ieta
    
    v = zeros(1,length(t));
    w = zeros(1,length(t));

    v(1) = 0;
    w(1) = 0;

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
    
    veta = v;
    weta = w;
    
    venv = zeros(1,Npieces);
    tenv = zeros(1,Npieces);
    for k=2:Npieces
        if eta(k)>=eta(k-1)
            [venv(k),tenv(k)] = max(veta(jon(k):jon(k+1)-1));
            tenv(k) = tenv(k)*dt+ton(k);
        else
            [venv(k),tenv(k)] = min(veta(jon(k):jon(k+1)-1));
            tenv(k) = tenv(k)*dt+ton(k);
        end
    end
   
    figure
    hold on
    plot(t,veta,'-b','linewidth',2);
    plot(t,Ieta,'-','Color',lightgray,'linewidth',1);
    plot(tenv,venv,'or','linewidth',1)
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    
    % Ordering of eta according to the increasing amplitude 
    
    etaord = eta;
    L = 0;
    for j=1:Npieces;
        L = L+1;
        for k=L:Npieces
            if etaord(k)<etaord(j)
                aux = etaord(j);
                etaord(j) = etaord(k);
                etaord(k) = aux;
            end
        end 
    end
    
    Ieta = zeros(1,length(t));      
    for l=1:length(ton)-1        
        Ieta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = eta(l);
    end
    Ieta(end) = Ieta(end-1);
    
   
    Ietaord = zeros(1,length(t));      
    for l=1:length(ton)-1        
        Ietaord(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = etaord(l);
    end
    Ietaord(end) = Ietaord(end-1);
    
    % Computation of the solution in response to Ietaord
    
    v = zeros(1,length(t));
    w = zeros(1,length(t));

    v(1) = 0;
    w(1) = 0;

    for j=1:length(t)-1
        k1v = -gL*pwlv(v(j),gL,gc,alpha)-g1*w(j)+D*Ietaord(j);
        k1w = (v(j)-w(j))/tau1;
        av = v(j)+k1v*dt;
        aw = w(j)+k1w*dt;
        k2v = -gL*pwlv(av,gL,gc,alpha)-g1*aw+D*Ietaord(j+1);
        k2w = (av-aw)/tau1;
        v(j+1) = v(j)+(k1v+k2v)*dt/2;
        w(j+1) = w(j)+(k1w+k2w)*dt/2;
    end
    
    vord = v;
    word = w;
    
    vordenv = zeros(1,Npieces);
    tordenv = zeros(1,Npieces);
    for k=2:Npieces
        if etaord(k)>=etaord(k-1)
            [vordenv(k),tordenv(k)] = max(vord(jon(k):jon(k+1)-1));
            tordenv(k) = tordenv(k)*dt+ton(k);
        else
            [vordenv(k),tordenv(k)] = min(vord(jon(k):jon(k+1)-1));
            tordenv(k) = tordenv(k)*dt+ton(k);
        end
    end
   
    figure
    hold on
    plot(t,vord,'-b','linewidth',2);
    plot(t,Ietaord,'-','Color',lightgray,'linewidth',1);
    plot(tordenv,vordenv,'or','linewidth',1);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    
    SEQ = 0;
            % 1: One input sequence (and all permutations) for each
            %    parameter set
           
            
    if SEQ == 1                 
    
        %   Sequence of random permutations of the order of constant pieces

        Ntrials = 1000;



        vpermsum = zeros(1,length(t));

        vpermenvord_vec = zeros(Ntrials,Npieces);
        tpermenvord_vec = zeros(Ntrials,Npieces);
        
        for l=1:Ntrials

            l 

            P = randperm(length(etaord),length(etaord));

            etaperm = etaord(P);
            Ietaperm = zeros(1,length(t));    
            for k=1:length(ton)-1        
                Ietaperm(floor(ton(k)/dt):floor((ton(k)+Tdur)/dt)) = etaperm(k);
            end
            Ietaperm(end) = Ietaperm(end-1);

            % Computation of the solution in response to Ietaperm

            v = zeros(1,length(t));
            w = zeros(1,length(t));

            v(1) = 0;
            w(1) = 0;
           
            for j=1:length(t)-1
                k1v = -gL*pwlv(v(j),gL,gc,alpha)-g1*w(j)+D*Ietaperm(j);
                k1w = (v(j)-w(j))/tau1;
                av = v(j)+k1v*dt;
                aw = w(j)+k1w*dt;
                k2v = -gL*pwlv(av,gL,gc,alpha)-g1*aw+D*Ietaperm(j+1);
                k2w = (av-aw)/tau1;
                v(j+1) = v(j)+(k1v+k2v)*dt/2;
                w(j+1) = w(j)+(k1w+k2w)*dt/2;
            end

            vperm = v;
            wperm = w;

            vpermenv = zeros(1,Npieces);
            tpermenv = zeros(1,Npieces);
            for k=2:Npieces
                if etaperm(k)>=etaperm(k-1)
                    [vpermenv(k),tpermenv(k)] = max(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                else
                    [vpermenv(k),tpermenv(k)] = min(vperm(jon(k):jon(k+1)-1));
                    tpermenv(k) = tpermenv(k)*dt+ton(k);
                end
            end

            if g1 == 0            
                vpermenv(1:length(ton)-1) = vperm(ton(2:end)/dt);
                tpermenv = ton(2:end);
            end
            
            vpermenvord_vec(l,P) = vpermenv;

            vpermsum = vpermsum+vperm;
            
            vpermenvordmax = zeros(1,Npieces);
            vpermenvordmin = zeros(1,Npieces);
            vpermenvordmean = zeros(1,Npieces);
            vpermenvordvar = zeros(1,Npieces);
            for k=1:Npieces
                vpermenvordmax(k) = max(vpermenvord_vec(:,k));
                vpermenvordmin(k) = min(vpermenvord_vec(:,k));
                vpermenvordmean(k) = mean(vpermenvord_vec(:,k));
                vpermenvordvar(k) = var(vpermenvord_vec(:,k));
            end
            
            
            
        end
        
        figure
        hold on
        plot(-100,-100,'-b','linewidth',2);
        %plot(-100,-100,'-r','linewidth',2);
        plot(-100,-100,'-g','linewidth',2);
        for l=1:Ntrials
            plot(vpermenvord_vec(l,:),'o');
        end
        plot(vpermenvordmean,'-b','linewidth',3)
        plot(vordenv(2:end),'-g','linewidth',3) 
        axis([0 Npieces -15 15]);
        set(gca,'fontsize',24);
        xlabel('Piece #');
        ylabel('');
        legend('<P_{\eta}>','P_{\eta,step}')


        figure
        hold on
        plot(t,vperm,'-b','linewidth',2);
        plot(tpermenv,vpermenv,'o','Color',lightcoral,'linewidth',2)
        set(gca,'fontsize',24);
        xlabel('t  [ms]');
        ylabel('V  [mV]');
        legend('V','P_{\eta}');

        figure
        hold on
        plot(vpermenvordvar,'-b','linewidth',2)
        plot(vpermenvordvar/max(va),'-r','linewidth',2)
        set(gca,'fontsize',24);
        xlabel('Piece #');
        ylabel('var(P{\eta})');
        legend('Var(P_{\eta})','VarN(P_{\eta})');   

        Data = vpermenvordvar';
        Datanorm = vpermenvordvar'/max(va);
        
 
    end
    
    
    PWS = 1;
            % 1: PSD and Z for a PWC signal with random amplitudes
            % 2: PSD and Z for a PWC signal with non-equal, but non-random
            %    amplitudes
            %
            % PWC: Piecewise constant
            % PSD: Power spectrum density 
            % Z:   Impedance
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


        % Random input amplitudes

        eta = randn(1,Npieces);
        
        Ieta = zeros(1,length(t));    
        for l=1:length(ton)-1        
            Ieta(floor(ton(l)/dt):floor((ton(l)+Tdur)/dt)) = eta(l);
        end
        Ieta(end) = Ieta(end-1);


        % Computation of the solution in response to Ieta

        v = zeros(1,length(t));
        w = zeros(1,length(t));

        v(1) = Ieta(1)/(gL+g1);
        w(1) = Ieta(1)/(gL+g1);
        v(1) = v(1)*1.2;
        
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


        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,D*Ieta);
        ZAP_Psd = PsdManSmooth;

        [PSD, Freq, PsdManSmooth, freqbin] = powerspectrum(t,v);
        V_Psd = PsdManSmooth;

        Z = V_Psd./ZAP_Psd;

        figure
        hold on
        plot(-100,-100,'ob','linewidth',2);
        plot(-100,-100,'or','linewidth',2);
        plot(freqbin,Z,'ob','linewidth',1);
        plot(freqbin,V_Psd*max(Z)/max(V_Psd),'or','linewidth',1);
        plot(fres,0,'ok','linewidth',3);
        plot(fnat,0,'ok','linewidth',3);
        plot([fres fres],[0 20],'--k','linewidth',2);
        plot([fnat fnat],[0 20],'--k','linewidth',2);
        %plot(freqbin,Z,'-b','linewidth',1);
        axis([0 100,-0.01,6]);
        set(gca,'fontsize',20);
        xlabel('Freq.  [Hz]');
        ylabel('Z');
        
%         figure
%         hold on
%         plot(-100,-100,'ob','linewidth',2);
%         plot(-100,-100,'or','linewidth',2);
%         plot(freqbin,V_Psd,'ob','linewidth',1);
%         %plot(freqbin,Z,'-b','linewidth',1);
%         axis([0 100,-0.01,2]);
%         set(gca,'fontsize',20);
%         xlabel('Freq.  [Hz]');
%         ylabel('PSD');
            
   end
    
end

if CNT == 1

    % Time definitions

    Tmax = 100000;
    dt = 0.1;
    t = 0:dt:Tmax;

    % White (Gaussian) noise

    eta = randn(1,length(t));

    


    % Computation of the numerical solution with white noise

    vs = zeros(1,length(t));
    ws = zeros(1,length(t));

    for j=1:length(t)-1
        k1v = -gL*vs(j)-g1*ws(j);
        k1w = (vs(j)-ws(j))/tau1;
        av = vs(j)+k1v*dt;
        av = av+sqrt(2*D*dt)*eta(j);
        aw = ws(j)+k1w*dt;
        k2v = -gL*av-g1*aw;
        k2w = (av-aw)/tau1;
        vs(j+1) = vs(j)+(k1v+k2v)*dt/2;
        vs(j+1) = vs(j+1)+sqrt(2*D*dt)*eta(j);
        ws(j+1) = ws(j)+(k1w+k2w)*dt/2;
    end

    figure
    hold on
    plot(t,vs,'-b','linewidth',2);
    axis([0 Tmax -4 5]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V');

    [mean(vs) var(vs)]
    
    PWS = 0;
    if PWS == 1
    
        [PSD,Freq, PsdManSmooth,freqbin] = powerspectrum(t,vs);
        V_Psd = PsdManSmooth;

        figure
        hold on
        plot(freqbin,V_Psd,'ob','linewidth',1);
        plot(freqbin,smooth(V_Psd,'lowess',13),'-b','linewidth',2);
        axis([0 100,-0.01,0.1]);
        set(gca,'fontsize',20);
        xlabel('Freq.  [Hz]');
        ylabel('PSD');

        [PSDwn,Freq,PsdwnManSmooth,freqbin] = powerspectrum(t,eta);
        WN_Psd = PsdwnManSmooth;
        WN_Psd = WN_Psd*sqrt(2);

        figure
        hold on
        plot(freq,Zanl,'-b','linewidth',2);
        plot(freqbin,V_Psd./WN_Psd,'or','linewidth',1);
        %plot(freqbin,V_Psd,'or','linewidth',1);
        plot([0 max(freq)],[0 0],':');
        axis([0 freq(end) 0 max(Zanl)*1.2])    
        set(gca,'fontsize',24);
        xlabel('f  [Hz]');
        ylabel('Z');
        legend('Z_{AUT}','Z_{WN}')
        title('Impedance profile')
    end
end


