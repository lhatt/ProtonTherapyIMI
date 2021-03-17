function lhd=likelihood(a,theta)
% Calculate the likelihood function for data points theta given a

x = linspace(0,1,1000);
dx = x(2)-x(1);
k = length(a);
n = 0:k;
u = n'*x;
a=[1 a];
p = (a/sum(abs(a)))*sin(pi*u)+1;
p = p/sum(p)/dx;   % Probability distribution for theta, which is
                   % the position the particle lands on the circle
       
% Find the index of the probability distribution that corresponds to the
% data point theta
for i=1:length(theta)
    index(i)=find(2*pi*x==theta(i));
end
% Calculate the likelihood function
lhd=prod(p(index));