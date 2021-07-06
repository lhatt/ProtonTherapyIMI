% A dose profile is calculated for 3 regions representing skin/bone/tumour.
% For a given dose profile, we estimate the 3 diffusion coefficients 
% within each region.
% A particle filter algorithm is used to solve this inverse problem. The 
% mutation step applies the Metropolis-Hastings method. 

clear all;
N=200;%No. of particles
ndp=1000;%No. of data points observed
m=30;%No. of iterations
h=0.2; %vertical distance from gamma particle origin (x) to detector 
no_b=6;%no. of bins/compartments for angle ranges
global max_iterations
max_iterations = 5;

% Set to true to plot particle positions at every step:
plot_cloud = false;

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

%Calculating the weights of each probability distribution
for i=1:no_b
    wPT(i)=trapz(xds,PTs(i,:));
end

%Observe ndp datapoints
for i=1:ndp
    bin_dps(i)=randsample(1:no_b,1,true,wPT);%Choose a bin
    dps(i)=randsample(xds,1,true,PTs(bin_dps(i),:));%Observe a datapoint
    index_dps(i)=find(xds==dps(i));%Record the index
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
    plot_particles(as,a);
end

% Initialise the covariance matrix
CV0=eye(k)*(.01/k);

% Setup matrix to store average 
av_as = zeros(m,length(a));

% Vector to plot l2 distance of particle filter to truth
l2_p = zeros(m,1);

for j=1:m
    % Calculate the importance weights using the likelihood function and
    % the observations
    clear w;
    j
    for i=1:N
        % Dose profile corresponding to particle i
        [xs,Px]=DoseProfile(as(i,1),as(i,2),as(i,3));
        Pis=Prob_xd_bins(xs,Px,h,no_b); % P(xd) split into no_b bins by angle range
        wgt=0;
        % Use Log to avoid underflow
        for ii=1:ndp
            wgt=wgt+log(Pis(bin_dps(ii),index_dps(ii)));
        end
        w(i)=wgt+660*log(10);%multiplying the weight by large number
        % Weight of particle i
        w(i)=exp(w(i));
        
    end
    
    w=w/sum(w);% Normalise the weights
    % Resample n particles according to the weights
    as_w=datasample(as,N,'Weights',w);
    
    % Plot particle cloud
    if plot_cloud
        plot_particles(as_w,a);
    end
    %Calculating the expected value for a
    for i=1:N
        av_as(j,:)=av_as(j,:)+ w(i)*as(i,:);
    end    
    % Calculate the L2 error
    l2_p(j) = (av_as(j,:)-a)*(av_as(j,:)-a)'
    av_as(j,:)
   % Mutation step
    clear as_m;
    for i=1:N
        %Propose new step
        a2=[-1,-1,-1];
        while sum(a2<0)>0
            a2=mvnrnd(as_w(i,:),CV0,1);
        end
        a1=as_w(i,:); % past step
        
        % Find the distributions of the past and new state

        [xs1,Px1]=DoseProfile(a1(1),a1(2),a1(3));%dose profile past state
        [xs2,Px2]=DoseProfile(a2(1),a2(2),a2(3));%dose profile new state

        Pi1s=Prob_xd_bins(xs1,Px1,h,no_b);
        Pi2s=Prob_xd_bins(xs2,Px2,h,no_b);
        
        % Calculate the likelihoods of the past and new state
        lh1=0;lh2=0;
        for ii=1:ndp
            lh1=lh1+log(Pi1s(bin_dps(ii),index_dps(ii)));
            lh2=lh2+log(Pi2s(bin_dps(ii),index_dps(ii)));
        end
        lh1=lh1+660*log(10);
        lh2=lh2+660*log(10);
        lh1=exp(lh1);
        lh2=exp(lh2);
        % Calculate the weights for the past and new state
        w1=lh1*exp(-(1/(2*stda^2))*sum((a1-mean_a).^2));
        w2=lh2*exp(-(1/(2*stda^2))*sum((a2-mean_a).^2));
        alpha=min([1, w2/w1]);%acceptance probability
        as_m(i,:)=datasample([a1;a2],1,'Weights',[1-alpha;alpha]);% select new step
    end      
    as=as_m;% adopt the 'mutated' particles
   %Observe ndp datapoints
    for i=1:ndp
        bin_dps(i)=randsample(1:no_b,1,true,wPT);
        dps(i)=randsample(xds,1,true,PTs(bin_dps(i),:));
        index_dps(i)=find(xds==dps(i));
    end 
end
% Plot log(L2 error)
figure;
plot(1:m,log(l2_p));
xlabel 'iteration';
ylabel 'log(L2 error)'
grid on;