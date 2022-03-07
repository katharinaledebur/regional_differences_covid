function [x,CI,resnorm,residual,exitflag,output,RSS] = run_solver(Bezirk,beta,t_delay,training_obs,x0,t_cal)
    
    % solve for variables to get effect sizes (levenberg-marquardt)
    RSS = zeros(2);
    options = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',1000, 'FunctionTolerance',0.001);% 
    options.Algorithm = 'levenberg-marquardt';
    [x,resnorm,residual,exitflag,output,~,jacobian] = lsqnonlin(@nestedfun,x0,[],[],options);
    CI = nlparci(x,residual,'jacobian',jacobian);
    
    
    function y = nestedfun(x0)
        
        N = numel(Bezirk);
        T = size(Bezirk(1).Cases,1);
        c_bez_mod = zeros(N,T,4); 
        c_bez_mod_0 = zeros(N,T,4); 
        c_bez_dat = zeros(N,T,4); 
        i_bez_dat = zeros(N,T,4);
        
        for i=1:N
            for a=1:4
                c_bez_dat(i,:,a) = (Bezirk(i).Cases(:,a))./Bezirk(i).N(a);
            end
        end
        for a=1:4
            i_bez_dat(:,:,a) = [zeros(116,1) diff(squeeze(c_bez_dat(:,:,a))')'];
        end

        % calculate null- and augmented-model for every district 
        for i=1:N    
            b = Bezirk(i).BL;
            S = zeros(T,4);
            I = zeros(T,4);
            R = zeros(T,4);
            alpha = Bezirk(i).alpha_0;
            t0 = find(sum(Bezirk(i).Cases,2),1);
            I(t0,:) = (Bezirk(i).Cases(t0,:))./Bezirk(i).N; S(t0,:) = 1-I(t0,:);
            cm = Bezirk(i).CM;
            for a=1:4
                c_bez_mod(i,t0,a) = c_bez_dat(i,t0,a);
                c_bez_mod_0(i,t0,a) = c_bez_mod(i,t0,a);
            end
            % null-model
            for t=t0+1:T
                for a=1:4
                    sia = S(t-1,a)*sum(cm(a,:).*I(t-1,:));
                    S(t,a) = S(t-1,a) - alpha(t-1,a)*sia;
                    I(t,a) = I(t-1,a) + alpha(t-1,a)*sia - beta*I(t-1,a);
                    R(t,a) = R(t-1,a) + beta*I(t-1,a);
                    c_bez_mod_0(i,t,a) = c_bez_mod_0(i,t-1,a) + alpha(t-1,a)*sia;
                end
            end
            % augmented model
            alpha_bez = zeros(T,4);
             for t=t0+1:T
                tm = max(t-t_delay,1); 
                for a=1:4
                    al = alpha(t-1,a) * (1 + x0(1)*Bezirk(i).Temp(tm)) * (1 + x0(2)*Bezirk(i).Cloud(tm)) * (1 + x0(3)*Bezirk(i).Humid(tm)) * (1 + x0(4)*Bezirk(i).Prec(tm)) * (1 + x0(5)*Bezirk(i).Wind(tm)); % Metereological factors
                    if a==1
                       al = al * (1 + x0(6)*Bezirk(i).Measures(1,tm)); % School < 20y
                    else
                        al = al * (1 + x0(7)*Bezirk(i).Measures(1,tm)); % School > 20y
                    end
                    al = al * (1 + x0(8)*Bezirk(i).Measures(2,tm)); % Gastronomy
                    al = al * (1 + x0(9)*Bezirk(i).Measures(3,tm)); % Healthcare
                    al = al * (1 + x0(10)*Bezirk(i).Measures(4,tm)); % Mass Events
                    al = al * (1 + x0(11)*Bezirk(i).Rg(tm));
                    alpha_bez(t,a) = al;
                    sia = S(t-1,a)*sum(cm(a,:).*I(t-1,:));
                    S(t,a) = S(t-1,a) - al*sia;
                    I(t,a) = I(t-1,a) + al*sia - beta*I(t-1,a);
                    R(t,a) = R(t-1,a) + beta*I(t-1,a);
                    c_bez_mod(i,t,a) = c_bez_mod(i,t-1,a) + al*sia;
                end
            end
        end
        i_bez_dat = zeros(N,sum(t_cal),4);
        for a=1:4
            i_bez_dat(:,:,a) = [zeros(116,1) diff(squeeze(c_bez_dat(:,t_cal,a))')'];
        end
        i_bez_mod = zeros(N,sum(t_cal),4);
        for a=1:4
            i_bez_mod(:,:,a) = [zeros(116,1) diff(squeeze(c_bez_mod(:,t_cal,a))')'];
        end
        i_bez_mod_0 = zeros(N,sum(t_cal),4);
        for a=1:4
            i_bez_mod_0(:,:,a) = [zeros(116,1) diff(squeeze(c_bez_mod_0(:,t_cal,a))')'];
        end
        % calculate residual sum of squares
        y = squeeze(reshape(i_bez_mod-i_bez_dat,1,1,[]));
        y0 = squeeze(reshape(i_bez_mod_0-i_bez_dat,1,1,[]));
        if sum(training_obs) == numel(training_obs)
            RSS(1,2) = sum(y0.^2);
            RSS(2,2) = sum(y.^2);
            RSS(1,1) = NaN; RSS(2,1) = NaN;
        else
            RSS(1,1) = sum(y0(training_obs==0).^2);
            RSS(1,2) = sum(y0(training_obs==1).^2);
            RSS(2,1) = sum(y(training_obs==0).^2);
            RSS(2,2) = sum(y(training_obs==1).^2);
        end
        
        y = y(training_obs>0);
    end
end