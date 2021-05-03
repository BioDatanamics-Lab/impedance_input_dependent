SEQ = 1;
            % 1: one input sequence (and all permutations) for each
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
        
        Data = vpermenvordvar';
         Datanorm = vpermenvordvar'/max(va);
            
            
            
            