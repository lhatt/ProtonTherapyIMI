function u = u_poly(x)

%u = x(1)*(1-x(1))*x(2)*(1-x(2));

u = exp(-1000*(x(1)^2+(x(2)-0.5)^2));