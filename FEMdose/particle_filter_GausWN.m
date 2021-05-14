% A dose profile is calculated for 3 regions representing skin/bone/tumour.
% For a given dose profile, we estimate the 3 diffusion coefficients 
% within each region.
% A particle filter algorithm is used to solve this inverse problem. The 
% mutation step applies the Metropolis-Hastings method. 

clear all;
N=100;%No. of particles
ndp=1;%No. of data points observed- one real dose profile
m=50;%No. of iterations
           
% Set to true to plot particl positions at every step:
plot_cloud = false;

% Find the true statistic
%Divide the domain into 3 chunks of tissue
% diffusion over region 1 (skin)
d1 = 1;
% diffusion over region 2 (bone)
d2 = 0.1;
% diffusion over region 3 (tumour)
d3 = 0.5;
a=[d1,d2,d3];% true statistic 
k=length(a);% length of true statistic
mean_a=mean(a);
stda=std(a);
%Observation data- dose profile
DP=Gdetector(d1,d2,d3);

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
sd=(2.38^2)/k;

% Setup matrix to store average 
av_as = zeros(m,length(a));

% Vector to plot l2 distance of particle filter to truth
l2_p = zeros(m,1);

for j=1:m
    % Calculate the importance weights using the likelihood function and
    % the data points
    clear w;
    for i=1:N
        % Data corresponding to particle i
        DPa=Gdetector(as(i,1),as(i,2),as(i,3));
        % Weight of particle i
        w(i)=likelihoodI(DPa',DP');
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
        av_as(j,:) =av_as(j,:)+ w(i)*as(i,:);
    end    
    % Calculate the L2 error
    l2_p(j) = (av_as(j,:)-a)*(av_as(j,:)-a)';
    
   % Mutation step
    clear as_m;
    for i=1:N
        %Propose new step
        a2=mvnrnd(as_w(i,:),CV0,1);
        a1=as_w(i,:); % past step
        % Calculate the likelihood of the past and new state
        DP1=Gdetector(a1(1),a1(2),a1(3));
        DP2=Gdetector(a2(1),a2(2),a2(3));
        % Calculate the weights for the past and new state
        w1=likelihoodI(DP1',DP')*exp(-(1/(2*stda^2))*sum((a1-mean_a).^2));
        w2=likelihoodI(DP2',DP')*exp(-(1/(2*stda^2))*sum((a2-mean_a).^2));
        alpha=min([1, w2/w1]);%acceptance probability
        as_m(i,:)=datasample([a1;a2],1,'Weights',[1-alpha;alpha]);% select new step
        
    end      
    as=as_m;
   
   
end
% Plot log(L2 error)
figure;
plot(1:m,log(l2_p));
xlabel 'iteration';
ylabel 'log(l2 error)'
grid on;