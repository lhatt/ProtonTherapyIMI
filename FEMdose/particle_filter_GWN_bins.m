% A dose profile is calculated for 3 regions representing skin/bone/tumour.
% For a given dose profile, we estimate the 3 diffusion coefficients 
% within each region.
% A particle filter algorithm is used to solve this inverse problem. The 
% mutation step applies the Metropolis-Hastings method, with the option of
% applying the Adaptive Metropolis technique.

clear all;
N = 200;%No. of particles
ndp = 1000;%No. of data points observed
m = 20;%No. of iterations
h = 0.2; %vertical distance from gamma particle origin (x) to detector 
no_b=6;%no. of bins/compartments for angle ranges
global max_iterations
max_iterations = 4;

id = datestr(now,'yymmddHHMMSS');
dir = 'output/';

filename_diary = strcat(id,'Diary.txt');

timerVal = tic;

diary(filename_diary)
diary on

% Set to true to plot particle positions at every step:
plot_cloud = true;

% Set to true if using adaptive covariance matrix for mutation step
AM = true;

% Find the true statistic
%Divide the domain into 3 chunks
% diffusion over region 1 (bone)
d1 = 1;
% diffusion over region 2 (tissue)
d2 = 0.1;
% diffusion over region 3 (tumour)
d3 = 0.5;
a=[d1,d2,d3];% true statistic: 3 diffusion coefficients
k=length(a);% length of true statistic
mean_a=mean(a);% mean of true statistic
stda=std(a);% s.d. of true statistic
%Observation data- True dose profile/prob. distn for x
[xs,PxT]=DoseProfile(d1,d2,d3);
%Calculate probability distribution P(xd)- split into no_b bins
%Each bin represents a different angle range
PTs=Prob_xd_bins(xs,PxT,h,no_b);
xds=xs(1:end-1);

% Store all data points
tot_dat = ndp*m;
dat = zeros(1,tot_dat);
bin_dat = zeros(1,tot_dat);
index_dat = zeros(1,tot_dat);

% Count acceptance rate for metropolis-hastings at each step.
accept_MH = zeros(1,m);

%Calculating the weights of each probability distribution
for i=1:no_b
    wPT(i)=trapz(xds,PTs(i,:));
end

fig = figure('Visible','off');
for i=1:6
    subplot(3,2,i)
    plot(xds,PTs(i,:))
end

saveas(fig,strcat(dir,id,'Truth'),'epsc');


%Observe ndp datapoints
for i=1:ndp
    bin_dat(i)=randsample(1:no_b,1,true,wPT);%Choose a bin
    index_dat(i)=randsample(1:length(xds),1,true,PTs(bin_dat(i),:));%Record the index
    dat(i)=xds(index_dat(i)); %Observe a datapoint
end
% Sampling N particles from Gaussian White noise Prior
clear as;
for i=1:N
    clear part;
    for j=1:k
        part(j)=-1;
        % Exclude negative values
        while part(j)<0
            part(j)=stda*randn + mean_a;
        end
    end
    as(i,:)=part;
end


%Plot the particles in state space
if plot_cloud
    fig = figure('Visible','off');
    plot_particles(as,a);
    saveas(fig,strcat(dir,id,'InitParts'),'epsc');
end

% Initialise the covariance matrix
CV0=eye(k)*(.01/k);
% Initialise adaptive covariance matrix
if AM
    for i=1:N
        CV{i}=CV0;
        CV_1{i}=as(i,:)'*as(i,:);
    end
    beta=0.05;
    sd=(2.38^2)/k;
    new_states=as';
end

% Setup matrix to store average 
av_as = zeros(m,length(a));

% Vector to plot l2 distance of particle filter to truth
l2_p = zeros(m,1);

% To store weight of particles after each iteration. A basic measure of how
% well the particles are predicting the outcomes.
log_wt = zeros(m,1);

% Vector to store the MSE
mse = zeros(m,1);

% Initialise dose profile vector
dp = zeros(N,length(xs))-2; % Will be negative if not initialised.

for j=1:m
    % Calculate the importance weights using the likelihood function and
    % the observations
    
    j
    
    % Compute log_weight of truth as a normalising constant
    wgt_normalise = 0;
    for ii=((j-1)*ndp+1):(j*ndp)
        wgt_normalise=wgt_normalise+log(PTs(bin_dat(ii),index_dat(ii)));
    end

    % Plot initial distributions.  
    fig = figure('Visible','off');
    hold on;
    
    for i=1:N
        % Dose profile corresponding to particle i
         if dp(i,1) < -1
             [xs,Px]=DoseProfile(as(i,1),as(i,2),as(i,3));
             dp(i,:)=Px;
         else
             [xs,Px]=DoseProfile(as(i,1),as(i,2),as(i,3));             
         end
        if dp(i,:)~=Px
            disp(as(i,:))
            dp(i,:) = Px
        end
        Pis=Prob_xd_bins(xs,dp(i,:),h,no_b); % P(xd) split into no_b bins by angle range

        wgt=0;
        % Use Log to avoid underflow
        for ii=((j-1)*ndp+1):(j*ndp)
            wgt=wgt+log(Pis(bin_dat(ii),index_dat(ii)));
        end
        w(i)=wgt-wgt_normalise; %multiplying the weight by large number
        % Weight of particle i
        w(i)=exp(w(i));
        
        %Plot initial profiles
        for kk=1:6
            %fig2 = 
            subplot(3,2,kk);
            hold on;
            plot(xds,Pis(kk,:),'Color',[1, 0, 0, 0.4],'LineWidth',1)
        end
    end
    
    %Plot truth
    for kk=1:6
        %fig2 = 
        subplot(3,2,kk);
        hold on;
        plot(xds,PTs(kk,:),'Color',[0, 0.4, 1, 1],'LineWidth',2)
    end
    
    filepath = char(strcat(dir,id,'Pred',string(j)));
    
    saveas(fig,filepath,'epsc');
    hold off;
    
    fig = figure('Visible','off');
    hold on;
    for kk = 1:N
        plot(xs,dp(kk,:),'Color',[0.4, 0, 0.2, 0.4],'LineWidth',1)
    end
    plot(xs,PxT,'Color',[0, 1, 0, 1],'LineWidth',2)

    filepath = char(strcat(dir,id,'Prof',string(j)));
    
    saveas(fig,filepath,'epsc');
    hold off;

    log_wt(j) = log(sum(w)/N)+wgt_normalise;
    % display(w)
    w=w/sum(w);% Normalise the weights
    
    % Resample n particles according to the weights
    perm = datasample(1:N,N,'Weights',w); % Setup permutation
    as_w=as(perm,:);
    dp = dp(perm,:);
    % Relabel the CV matrices
    CV_temp{N} = CV0; % Initialise matrix
    for i = 1:N
        CV_temp{i} = CV_1{perm(i)};
    end
    CV_1 = CV_temp;
    new_states = new_states(:,perm);        
        
    % Plot particle cloud
    if plot_cloud
        fig = figure('Visible','off');
        plot_particles(as_w,a);
        
        filepath = char(strcat(dir,id,'Part',string(j)));
        saveas(fig,char(filepath),'epsc');
    end
    %Calculating the expected value for a
    for i=1:N
        av_as(j,:)=av_as(j,:)+ w(i)*as(i,:);
    end    
    % Calculate the L2 error
    l2_p(j) = (av_as(j,:)-a)*(av_as(j,:)-a)';
    disp('L2 error in mean')
    disp(l2_p(j))
    disp('Average values of parameters')
    disp(av_as(j,:))

    % Calculate the MSE
    mse(j) = trace((as_w-a)'*(as_w-a))/N;
    disp('Mean Squared Error of Particles')
    disp(mse(j))
    
    % Mutation step
    clear as_m;
    for i=1:N
        % Compute the adaptive covariance matrix
        
        if AM
            CV_1{i}=as_w(i,:)'*as_w(i,:)+CV_1{i};
            new_states(:,i)=as_w(i,:)'+new_states(:,i);
            mn_st=new_states(:,i)/(j+1);
            CV{i}=(sd/j)*( CV_1{i}-(j+1)*mn_st*mn_st')+sd*.1*eye(k);
        end    
       
        %Propose new step
        a2=[-1,-1,-1];
        while sum(a2<0)>0
            if j>2*k && AM
                % Adaptive Metropolis step
                % Check to see if CV{i} is positive definite. If not,
                % default to CV0
                try
                    a2=(1-beta)*mvnrnd(as_w(i,:),CV{i},1)+...
                        beta*mvnrnd(as_w(i,:),CV0,1);
                catch ME
                    if (strcmp(ME.identifier,'stats:mvnrnd:BadCovariance2DSymPos'))
                        a2 = mvnrnd(as_w(i,:),CV0,1);
                        disp('Warning, non-positive definite covariance matrix')
                        i
                        CV{i}
                    else
                        rethrow(ME)
                    end
                end
            else
                a2=mvnrnd(as_w(i,:),CV0,1);
            end
            % a2=mvnrnd(as_w(i,:),CV0,1);
        end
        a1=as_w(i,:); % past step
        
        % Find the distributions of the past and new state

        [xs2,Px2]=DoseProfile(a2(1),a2(2),a2(3));%dose profile new state

        Pi1s=Prob_xd_bins(xs2,dp(i,:),h,no_b);
        Pi2s=Prob_xd_bins(xs2,Px2,h,no_b);
        
        % Calculate the likelihoods of the past and new state
        lh1=0;lh2=0;
        for ii=1:(j*ndp)
            lh1=lh1+log(Pi1s(bin_dat(ii),index_dat(ii)));
            lh2=lh2+log(Pi2s(bin_dat(ii),index_dat(ii)));
        end
        %lh1=lh1+660*log(10);
        %lh2=lh2+660*log(10);
        %lh1=exp(lh1);
        lh2=exp(lh2-lh1);
        lh1=1;
        % Calculate the weights for the past and new state
        w1=lh1*exp(-(1/(2*stda^2))*sum((a1-mean_a).^2));
        w2=lh2*exp(-(1/(2*stda^2))*sum((a2-mean_a).^2));
        alpha=min([1, w2/w1]);%acceptance probability
        % display(alpha);
        as_m(i,:)=datasample([a1;a2],1,'Weights',[1-alpha;alpha]);% select new step
        
        if as_m(i,:) == a2
            accept_MH(j) = accept_MH(j) + 1; % Count if we accepted
            dp(i,:) = Px2; % Update stored profile
        end
        
    end      
    as=as_m;% adopt the 'mutated' particles
    
    %Observe ndp datapoints
    for i=(j*ndp+1):((j+1)*ndp)
        bin_dat(i)=randsample(1:no_b,1,true,wPT);
        index_dat(i)=randsample(1:length(xds),1,true,PTs(bin_dat(i),:));%Record the index
        dat(i)=xds(index_dat(i)); %Observe a datapoint
    end 
end
% Plot log(L2 error)
fig = figure();
plot(1:m,log(l2_p));
xlabel 'iteration';
ylabel 'log(L2 error)'
grid on;

filepath = char(strcat(dir,id,'Error'));

saveas(fig,filepath,'epsc');

fig = figure();
plot(1:m,accept_MH/N);
xlabel 'iteration';
ylabel 'Proportion of accepted MH proposals'
grid on;

filepath = char(strcat(dir,id,'MHAccept'));

saveas(fig,filepath,'epsc');

fig = figure();
plot(1:m,log_wt);
xlabel 'iteration';
ylabel 'Weight of particles after observing data'
grid on;

filepath = char(strcat(dir,id,'weight'));

saveas(fig,filepath,'epsc');

fig = figure();
plot(1:m,mse);
xlabel 'iteration';
ylabel 'MSE of particles'
grid on;

filepath = char(strcat(dir,id,'MSE'));

saveas(fig,filepath,'epsc');

toc(timerVal)

diary off
movefile(filename_diary,dir)

