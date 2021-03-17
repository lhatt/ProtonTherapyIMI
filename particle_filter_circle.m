% A source that releases particles is at the center of a circular detector,
% which records the location of the particles.
% A particle filter algorithm is used to solve this inverse problem. The 
% mutation step applies the Metropolis-Hastings method. As well, there is
% an option to use the Adaptive Metropolis technique for an adaptive normal 
% distribution in the mutation step.

clear all;
N=1000;%No. of particles
k=3;%Length of summary statistic
ndp=500;%No. of data points observed per iteration
m=5;%No. of total iterations
AM_start=2;%When to introduce adaptive normal distribution in mutation step
           %Set AM_start>m for a fixed covariance matrix 

% Find the true statistic
a=sum_stat(k);
max_a=linspace(1,k+2,k+1).^(-1);
max_a=max_a(2:end);%Bounds for a vector

% Calculate the corresponding probability distribution on the circle
x = linspace(0,1,1000);
dx = x(2)-x(1);
n = 0:k;
u = n'*x;
ah=[1 a];
p = (ah/sum(abs(ah)))*sin(pi*u)+1;
p = p/sum(p)/dx;%Probability distribution for true stat.
theta = 2*pi*randsample(x,ndp,true,p*dx);%Observing ndp data points

% Sampling N particles from distribution for a
clear as;
for i=1:N
    as(i,:)=sum_stat(k);
end
%Plot the particles in state space
plot_particles(as,a);

% Initialise the covariance matrix
for i=1:N
    CV{i}=eye(k)*.01;
end
sum_states=as';


% Introduce ndp data points per cycle
for j=1:m
    % Calculate the importance weights using the likelihood function and
    % the ndp data points
    clear w;
    for i=1:N
        w(i)=likelihood(as(i,:),theta);
    end
    w=w/sum(w);% Normalise the weights
    
    % Resample n particles according to the weights
    as_w=datasample(as,N,'Weights',w);
    plot_particles(as_w,a);
     
   % Mutation step
    clear as_m;
    acc_count=0;
    for i=1:N
        sum_states(:,i)=as_w(i,:)'+sum_states(:,i);
        % Calculate the mean state
        mn_st=sum_states(:,i)/(j+1);
        % Calculate the covariance matrix using Adaptive Metropolis alg.
        if j>=AM_start
             CV{i}=((j-1)/j)*CV{i}+(1/j)*(as_w(i,:)'-mn_st)*...
                 (as_w(i,:)'-mn_st)';
        end
        % Sample new state from normal distribution with covariance matrix
        % CV, its mean is the previous particle
        a2=max_a+.1;%Ensure a is bounded 
        while sum(abs(a2)>max_a)>0
            a2=mvnrnd(as_w(i,:),CV{i},1);
        end
        
        a1=as_w(i,:);
        % Calculate the likelihood of the past and new states
        w1=likelihood(a1,theta);
        w2=likelihood(a2,theta);
        alpha=min([1, w2/w1]);%acceptance probability
        as_m(i,:)=datasample([a1;a2],1,'Weights',[1-alpha;alpha]);
        if as_m(i,:)==a2
            acc_count=acc_count+1;
        end
    end      
    acc_rate(j)=acc_count/N;%acceptance rate
    as=as_m;
   
    %Make ndp new observations for theta
    x = linspace(0,1,1000);
    theta = 2*pi*randsample(x,ndp,true,p*dx);

end

% Calculate cluster characteristics
mean_acc_rate=mean(acc_rate);%Mean acceptance rate
centroid=mean(as);%Centroid of the cluster
dist=norm(centroid-a);%Distance between centroid and true stat.
for i=1:N
    sci(i)=norm(centroid-as(i,:));
end
sc=sum(sci)/N;%Size of the cluster
varW= (sum(w.^2))^(-1);%Variance of the weights or ESS

figure;
hold on;
for i=1:N
    lh = plot(x,ThetaDensity(x,as(i,:)));
    lh.Color = [1,0,0,0.5];
end
plot(x,ThetaDensity(x,a),'b','LineWidth',4);
hold off;


