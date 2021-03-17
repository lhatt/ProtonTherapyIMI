function p = ThetaDensity(x,a)
%THETADENSITY return the density on the circle of a given series, a; 
%             x is the "circle"

dx = x(2)-x(1);
k = length(a);
n = 0:k;
u = n'*x;
ah=[1 a];
p = (ah/sum(abs(ah)))*sin(pi*u)+1;
p = p/sum(p)/dx;%Probability distribution for true stat.
end

