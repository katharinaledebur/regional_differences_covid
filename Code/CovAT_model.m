clear all
close all

%% IMPORT CASE DATA, POPULATION AND CONTACT MATRICES

C = readtable('cum_cases_district_age.csv'); % for everyday, district and age group the cumulative covid cases (488x4 - from 2020-02-26 to 2021-06-27)
y = readtable('population_age.csv','HeaderLines',2); % population for every district and age group

% create struct to collect all of the data inouts for every district, date
% and age group
Bezirk  = struct;
[uID, I, J] = unique(C.BezirkID);
tau = min(C.date):max(C.date);

c_at = zeros(numel(tau),9,4); n_at = zeros(9,4);
for i=1:numel(uID)
    ind = J==i;
    Bezirk(i).ID = uID(i);
    Bezirk(i).Cases = zeros(numel(tau),4);
    x = C(ind,:);
    [~,k] = ismember(tau,x.date);
    Bezirk(i).Cases(find(k),1) = x.cumulative_sum_0;
    Bezirk(i).Cases(find(k),2) = x.cumulative_sum_1;
    Bezirk(i).Cases(find(k),3) = x.cumulative_sum_2;
    Bezirk(i).Cases(find(k),4) = x.cumulative_sum_3;
    k = y{:,2} == uID(i);
    Bezirk(i).N = y{find(k),3:6};
    b = floor(uID(i)/100);
    Bezirk(i).BL = b;
    for a=1:4
        c_at(:,b,a) = c_at(:,b,a) + Bezirk(i).Cases(:,a);
        n_at(b,a) = n_at(b,a) + Bezirk(i).N(a);
    end
    cm = readtable(['contact_matrices/',num2str(uID(i)),'.csv'],'HeaderLines',1);
    Bezirk(i).CM = table2array(cm);
    for a = 1:4
        Bezirk(i).CM(a,:) = Bezirk(i).CM(a,:)/sum(Bezirk(i).CM(a,:));
    end
end
% calculation of the contact matrices for every federal state (for
% calculation of the federal state transmission rates)
CMb = struct;
for b=1:9
    ind = find([Bezirk.BL]==b);
    x = zeros(4); n = 0;
    for i=ind
        x = x + Bezirk(i).CM.*sum(Bezirk(i).N);
        n = n + sum(Bezirk(i).N);
    end
    CMb(b).x = x./n;
end

%% IMPORT WEATHER DATA
% format the input timeseries (select 50% percentile) and collect
% everything in Bezirk struct
t_ini = datetime(2020,7,1);
t_out = datetime(2021,5,15);

t_cal = ismember(tau,t_ini:t_out);

stds = [];
avg_weather = zeros(365,5);
x = table2cell(readtable('temp_mean_Bezirk.csv'));
y = zeros(numel(uID),numel(tau));
tx = datetime(x(3:end,1));
[tt,It] = ismember(tau,tx);
for i=1:numel(uID)
    ind1 = strcmp(x(2,:),'50%');
    ind2 = strcmp(x(1,:),num2str(uID(i)));
    m = find(ind1.*ind2);
    y(i,tt) = str2double(x(It(It>0)+2,m));
end
temp = y;

for t=1:numel(tau)
    for i=1:numel(uID)
        indb = [Bezirk.BL]==Bezirk(i).BL; indb(i) = 0;
        temp(i,t) = y(i,t)-mean(y(indb,t));
    end
end
for i=1:numel(uID)
    Bezirk(i).Temp = temp(i,:);
end
for t=1:365, avg_weather(t,1) = mean(str2double(x(t+2,ind1))); end
stds(1) = std(reshape(temp(:,t_cal),1,[]));

x = table2cell(readtable('cloud_mean_Bezirk.csv'));
y = zeros(numel(uID),numel(tau));
tx = datetime(x(3:end,1));
[tt,It] = ismember(tau,tx);
for i=1:numel(uID)
    ind1 = strcmp(x(2,:),'50%');
    ind2 = strcmp(x(1,:),num2str(uID(i)));
    m = find(ind1.*ind2);
    y(i,tt) = str2double(x(It(It>0)+2,m));
end
cloud = y;
for t=1:numel(tau)
    for i=1:numel(uID)
        indb = [Bezirk.BL]==Bezirk(i).BL; indb(i) = 0;
        cloud(i,t) = y(i,t)-mean(y(indb,t));
    end
end
for i=1:numel(uID)
    Bezirk(i).Cloud = cloud(i,:);
end
for t=1:365, avg_weather(t,2) = mean(str2double(x(t+2,ind1))); end
stds(2) = std(reshape(cloud(:,t_cal),1,[]));

x = table2cell(readtable('humid_mean_Bezirk.csv'));
y = zeros(numel(uID),numel(tau));
tx = datetime(x(3:end,1));
[tt,It] = ismember(tau,tx);
for i=1:numel(uID)
    ind1 = strcmp(x(2,:),'50%');
    ind2 = strcmp(x(1,:),num2str(uID(i)));
    m = find(ind1.*ind2);
    y(i,tt) = str2double(x(It(It>0)+2,m));
end
humid = y;
for t=1:numel(tau)
    for i=1:numel(uID)
        indb = [Bezirk.BL]==Bezirk(i).BL; indb(i) = 0;
        humid(i,t) = y(i,t)-mean(y(indb,t));
    end
end
for i=1:numel(uID)
    Bezirk(i).Humid = humid(i,:);
end
for t=1:365, avg_weather(t,3) = mean(str2double(x(t+2,ind1))); end
stds(3) = std(reshape(humid(:,t_cal),1,[]));

x = table2cell(readtable('prec_mean_Bezirk.csv'));
y = zeros(numel(uID),numel(tau));
tx = datetime(x(3:end,1));
[tt,It] = ismember(tau,tx);
for i=1:numel(uID)
    ind1 = strcmp(x(2,:),'50%');
    ind2 = strcmp(x(1,:),num2str(uID(i)));
    m = find(ind1.*ind2);
    y(i,tt) = str2double(x(It(It>0)+2,m));
end
prec = y;
for i=1:numel(uID)
    Bezirk(i).Temp0 = prec(i,:);
end
for t=1:numel(tau)
    for i=1:numel(uID)
        indb = [Bezirk.BL]==Bezirk(i).BL; indb(i) = 0;
        prec(i,t) = y(i,t)-mean(y(indb,t));
    end
end
prec(isnan(prec)) = 0;
for i=1:numel(uID)
    Bezirk(i).Prec = prec(i,:);
end
for t=1:365, avg_weather(t,4) = mean(str2double(x(t+2,ind1))); end
stds(4) = nanstd(reshape(prec(:,t_cal),1,[]));

x = table2cell(readtable('wind_mean_Bezirk.csv'));
y = zeros(numel(uID),numel(tau));
tx = datetime(x(3:end,1));
[tt,It] = ismember(tau,tx);
for i=1:numel(uID)
    ind1 = strcmp(x(2,:),'50%');
    ind2 = strcmp(x(1,:),num2str(uID(i)));
    m = find(ind1.*ind2);
    y(i,tt) = str2double(x(It(It>0)+2,m));
end
wind = y;
for t=1:numel(tau)
    for i=1:numel(uID)
        indb = [Bezirk.BL]==Bezirk(i).BL; indb(i) = 0;
        wind(i,t) = y(i,t)-mean(y(indb,t));
    end
end
for i=1:numel(uID)
    Bezirk(i).Wind = wind(i,:);
end
for t=1:365, avg_weather(t,5) = mean(str2double(x(t+2,ind1))); end
stds(5) = std(reshape(wind(:,t_cal),1,[]));

wlab = {'Temperature', 'Cloudiness', 'Humidity','Precipitation', 'Wind'};

%% IMPORT MEASURES
M = readtable('Edited_Massnahmen.csv');
% combination of the measures
mID = [2 8 10 18];
mlab = cell(size(mID));
M.measureID(M.measureID==4) = 2;
M.measureID(M.measureID==9) = 8;
M.measureID(M.measureID==12) = 10;
M.measureID(M.measureID==19) = 18;
M.measureID(M.measureID==20) = 18;
M.measureID(M.measureID==15) = 14;

ind = (M.Start==datetime(2020,2,1)).*(M.Ende==datetime(2021,4,5));
M(ind>0,:) = [];
for i=1:size(M,1)
    if M.measureID(i)==2 && M.Ende(i)>datetime(2020,12,23)
        M.Ende(i)=datetime(2020,12,23);
    end
end

ld2 = datetime(2020,11,17):datetime(2020,12,6);
[~, Ild] = ismember(ld2,tau);

for i=1:numel(uID)
    m = zeros(numel(mID),numel(tau));
    x = M(M.BezirkID==uID(i),:);
    
    for j=1:size(x,1)
        [k, I] = ismember(x.measureID(j),mID);
        if k
            tt = x.Start(j):x.Ende(j);
            [~, It] = ismember(tt,tau);
            m(I,nonzeros(It)) = 1;
            if isempty(mlab{I}), mlab(I) = x.measure(j); end
        end
    end
    m(1,Ild) = 0;
    Bezirk(i).Measures = m;
end

nm = zeros(numel(mID),1);
for i=1:numel(mID)
    m = zeros(numel(uID),sum(t_cal));
    for j=1:numel(uID)
        m(j,:) = Bezirk(j).Measures(i,t_cal);
    end
    stds(i+5) = std(m(:));
    nm(i) = sum(sum(diff(m')==1));
end

%% IMPORT MOBILITY (Radius of Gyration)

mob = load('rg_Bezirk.mat','rg');
rg = mob.rg ;
for i=1:numel(uID)
    for t=1:numel(tau)
        Bezirk(i).Rg(t) = rg(i,t)  ;
    end
end

stds(5+numel(mlab)+2) = std(reshape(log(rg(:,t_cal)),1,[])); 

%% MAIN MODEL + RESULT
% beta weights according to gamma distribution from: 
% Paul, S., Lorin, E., Estimation of COVID-19 recovery and decease periods in Canada using delay model. Scientific Reports. 2021, 11, 23763, https://doi.org/10.1038/s41598-021-02982-w
bw = load('beta_weights.mat');
beta_scan = zeros(28,11);
N = numel(Bezirk);
% for every beta (each recovery time) calculation of the effect sizes
for rec_time = 14:42
   
    % calculation of the transmission rate observed in all other districts of the same federal state
    for b=1:numel(Bezirk)
        alpha = zeros(numel(tau),4);
        beta = 1/rec_time;
        S = ones(numel(tau),4);
        I = zeros(numel(tau),4);
        R = zeros(numel(tau),4);
        bl = Bezirk(b).BL;
        C = squeeze(c_at(:,bl,:));
        C = C-Bezirk(b).Cases;
        C = C./(n_at(bl,:)-Bezirk(b).N);

        Cm = zeros(numel(tau),4);
        t0 = find(sum(C,2),1);
        I(t0,:) = (squeeze(c_at(t0,bl,:))'-Bezirk(b).Cases(t0,:))./(n_at(bl,:)-Bezirk(b).N); 
        S(t0,:) = 1-I(t0,:);
        cm = CMb(bl).x;
        
        for t=t0+1:numel(tau)
            for a=1:4
                sia = S(t-1,a)*sum(cm(a,:).*I(t-1,:));
                alpha(t-1,a) = (C(t,a)-C(t-1,a))/(sia);
                S(t,a) = S(t-1,a) - alpha(t-1,a)*sia;
                I(t,a) = I(t-1,a) + alpha(t-1,a)*sia - beta*I(t-1,a);
                R(t,a) = R(t-1,a) + beta*I(t-1,a);
                Cm(t,a) = Cm(t-1,a) + alpha(t-1,a)*sia;
            end
        end
        Bezirk(b).alpha_0 = alpha;
    end
    training_obs = ones(N,sum(t_cal),4);
    training_obs = squeeze(reshape(training_obs,1,1,[]));
    x0 = zeros(11,1);
    % calculation of the null + augmented model
    [Result.alphas,Result.CI,Result.resnorm,Result.residual,Result.exitflag,Result.output,Result.residual_0] = run_solver(Bezirk,beta,0,training_obs,x0,t_cal);
    
    vals = Result.alphas *100;
    for i=[1:5 11]
        vals(i,:) = vals(i,:).*stds(i); 
    end 
    beta_scan(rec_time-13,:) = vals;
   
end
% calculation of the weighted (according to gamma distribution) average of
% the effect sizes and the corresponding standard deviation
weighted_means = zeros(1,11);
for i = 1:11
    weighted_means(1,i) = sum(beta_scan(:,i).*bw.w')/sum(bw.w);
end

errorbars = zeros(1,11);
for j = 1:11     
    errorbars(1,j) = sqrt(sum((beta_scan(:,j) - weighted_means(1,j)).^2 .*(bw.w(:)./sum(bw.w))));
end
% PLOT SUMMARY OF EFFECT SIZES (Fig3)
vals = weighted_means;
labels = {'Temperature', 'Cloudiness', 'Humidity', 'Precipitation', 'Wind','Schools (<20y)','Schools (\geq20y)','Gastronomy','Healthcare','Mass events','Radius of \newlineGyration'};
figure;
set(gcf,'Position',[80    130    800    750]);
axes('Position',[0.12    0.24    0.35    0.7]);
errorbar(1:5,vals(1:5),errorbars(1:5),errorbars(1:5),'d','MarkerFaceColor',[0.2 0.4 0.8],'MarkerEdgeColor',[0.2 0.4 0.8],'MarkerSize',10,'CapSize',18,'LineWidth',2,'Color',[0.2 0.4 0.8]);
hold on
ylim([-50 40]);xlim([0.5 5.5]);
title('Weather','FontSize',24,'Color',[0.2 0.4 0.8]);
for i=1:4,  line([i i]+.5,[-80 40],'Color',[0.8 0.9 1],'LineWidth',3); end
set(gca,'FontSize',18,'XColor',[0.2 0.4 0.8],'YColor',[0.2 0.4 0.8],'XTick',1:5,'XTickLabel',labels(1:5));
plot([0 6],[0 0],'k:');
ylabel('Impact on transmission rate \alpha [%]','Color','k');

box on
xtickangle(90)

axes('Position',[0.48    0.24    0.35    0.7]);
errorbar(1:5,vals(6:10),errorbars(6:10),errorbars(6:10),'d','MarkerFaceColor',[0 0.4 0],'MarkerEdgeColor',[0 0.4 0],'MarkerSize',10,'CapSize',18,'LineWidth',2,'Color',[0 0.4 0]);
hold on
ylim([-50 40]);xlim([0.5 5.5]);
plot([0 6],[0 0],'k:');
for i=1:4
    line([i i]+.5,[-80 40],'Color',[0.9 1 0.8],'LineWidth',3);
end
title('Restrictions','FontSize',24,'Color',[0 0.4 0]);
set(gca,'YTickLabel','','FontSize',18,'XColor',[0 0.4 0],'YColor',[0 0.4 0],'XTick',1:5,'XTickLabel',labels(6:10));
box on
xtickangle(90)

axes('Position',[0.84    0.24    0.1    0.7]);
errorbar(1,vals(11),errorbars(11),errorbars(11),'d','MarkerFaceColor',[0.6 0.3 0.8],'MarkerEdgeColor',[0.6 0.3 0.8],'MarkerSize',10,'CapSize',18,'LineWidth',2,'Color',[0.6 0.3 0.8]);
hold on
ylim([-50 40]);xlim([0.5 1.5]);
plot([0 2],[0 0],'k:');
title('Mobility','FontSize',24,'Color',[0.6 0.3 0.8]);
set(gca,'YTickLabel','','FontSize',18,'XColor',[0.6 0.3 0.8],'YColor',[0.6 0.3 0.8],'XTick',1,'XTickLabel',labels(11));
box on
xtickangle(90)

%% CROSS-VALIDATED HYPERPARAMETER SEARCH
N_cv = 100;
CV = struct;
N = numel(Bezirk);
T = size(Bezirk(1).Cases,1);
        
rec_time_range = [4:35];
t_delay_range = [0];
% split data (time) into test and training data (20/80 %) and run model + calculate residual sum of squares
for rec_time = rec_time_range
    
    disp(num2str(rec_time));
    for b=1:numel(Bezirk)
        alpha = zeros(numel(tau),4);

        beta = 1/rec_time;
        S = ones(numel(tau),4);
        I = zeros(numel(tau),4);
        R = zeros(numel(tau),4);
        bl = Bezirk(b).BL;
        C = squeeze(c_at(:,bl,:));
        C = C-Bezirk(b).Cases;
        C = C./(n_at(bl,:)-Bezirk(b).N);

        Cm = zeros(numel(tau),4);
        t0 = find(sum(C,2),1);
        I(t0,:) = (squeeze(c_at(t0,bl,:))'-Bezirk(b).Cases(t0,:))./(n_at(bl,:)-Bezirk(b).N); 

        S(t0,:) = 1-I(t0,:);
        cm = CMb(bl).x;
        for t=t0+1:numel(tau)
            for a=1:4
                sia = S(t-1,a)*sum(cm(a,:).*I(t-1,:));
                alpha(t-1,a) = (C(t,a)-C(t-1,a))/(sia);
                S(t,a) = S(t-1,a) - alpha(t-1,a)*sia;
                I(t,a) = I(t-1,a) + alpha(t-1,a)*sia - beta*I(t-1,a);
                R(t,a) = R(t-1,a) + beta*I(t-1,a);
                Cm(t,a) = Cm(t-1,a) + alpha(t-1,a)*sia;
            end
        end
        Bezirk(b).alpha_0 = alpha;
    end
    for t_delay=t_delay_range
        for cv=1:N_cv
            training_sample = 0.8;
            training_obs = zeros(N,sum(t_cal),4);
            biweek_ind = floor((datenum(tau(t_cal))-min(datenum(tau(t_cal))))/28);
            [unique_week_ind,I,J] = unique(biweek_ind);
            for i=1:numel(Bezirk)
                for a=1:4
                    j = unique_week_ind(randperm(numel(unique_week_ind)));
                    j = j(1:floor(numel(j)*training_sample));
                    training_obs(i,:,a) = ismember(J,j);
                end
            end
            training_obs = squeeze(reshape(training_obs,1,1,[]));

            if cv > 1
                x0 = CV(1,rec_time,t_delay+1).alphas;
            else
                if rec_time == rec_time_range(1)
                    x0 = zeros(11,1);
                else
                    x0 = CV(1,rec_time_range(find(rec_time_range==rec_time)-1),t_delay+1).alphas;
                end
            end
            [CV(cv,rec_time,t_delay+1).alphas,CV(cv,rec_time,t_delay+1).CI,CV(cv,rec_time,t_delay+1).resnorm,CV(cv,rec_time,t_delay+1).residual,CV(cv,rec_time,t_delay+1).exitflag,CV(cv,rec_time,t_delay+1).output, CV(cv,rec_time,t_delay+1).residual_0] = run_solver(Bezirk,beta,t_delay,training_obs,x0,t_cal);
            CV(cv,rec_time,t_delay+1).training_obs = training_obs; 
            CV(cv,rec_time,t_delay+1).alphas';
        end
        
    end
end
rss_train = []; rss_0_train = []; 
rss_test = []; rss_0_test = []; 
d_rss_train = []; d_rss_test = [];
d_rss_train_std = []; d_rss_test_std = [];
wr1 = 1;
for rec_time = rec_time_range
    wr2 = 1;
    for t_delay=t_delay_range
        x = zeros(N_cv,1); y = zeros(N_cv,1);
        for cv=1:N_cv
            x(cv) = CV(cv,rec_time,t_delay+1).residual_0(2,1);
            y(cv) = CV(cv,rec_time,t_delay+1).residual_0(1,1);
        end
        rss_test(wr1,wr2) = mean(x);
        rss_0_test(wr1,wr2) = mean(y);
        d_rss_test(wr1,wr2) = mean((y-x)./y);
        d_rss_test_std(wr1,wr2) = std((y-x)./y);
        x = zeros(N_cv,1); y = zeros(N_cv,1);
        for cv=1:N_cv
            x(cv) = CV(cv,rec_time,t_delay+1).residual_0(2,2);
            y(cv) = CV(cv,rec_time,t_delay+1).residual_0(1,2);
        end
        rss_train(wr1,wr2) = mean(x);
        rss_0_train(wr1,wr2) = mean(y);
        d_rss_train(wr1,wr2) = mean((y-x)./y);
        d_rss_train_std(wr1,wr2) = std((y-x)./y);
        wr2 = wr2+1;
    end
    wr1 = wr1 + 1;
end

% PLOT HYPERPARAMETERSEARCH / CROSS-VALIDATION RESULTS (Fig2)
rec_time_range = [4:35];
col = parula(10);
figure;
errorbar(rec_time_range,100*d_rss_train,100*d_rss_train_std,'-d','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',col(3,:),'Color',col(3,:));
hold on
errorbar(rec_time_range,100*d_rss_test,100*d_rss_test_std,'-s','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',col(9,:),'Color','k');
ylabel('\Delta RSS [%]'); xlabel('recovery time [d], 1/\beta');
legend('String',{'training','test'},'Location','Best');