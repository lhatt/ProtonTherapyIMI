function P=Prob_xd_bins(xs,Prx,h,no_b)
% Calculate the probability distribution for xd (the position where the
% gamma particle lands on the detector).
% Split the distribution into different angle ranges (the angle which the
% gamma particle leaves the x axis).
dx=xs(2)-xs(1);
xds=xs(1:end-1);
dth=(pi/2)/(no_b/2);%size of the bins
thetas=-pi/2:dth:pi/2;%different angle ranges
for i=1:length(xds)
    xd=xds(i);
    %Compute P(x and xd) and then average over x to find P(xd)
    int=0;
    for j=1:length(xs)
        x=xs(j);
        
        Pxd_x=(1/(2*pi))*(atan((xd+dx-x)/h)-atan((xd-x)/h));%P(xd|x)
        Px=Prx(j);%P(x)
        theta(j)=atan((xd-x)/h);% angle given xd and x
        intj(j)=Pxd_x*Px*dx;% averaging over x
    end
    % Distribute each theta into the appropriate bin
    for j=1:no_b
        ind{j}=find((theta>=thetas(j)).*(theta<thetas(j+1)));
    end
    % Compute the integral for each bin
    for j=1:no_b
        P(j,i)=sum(intj(ind{j}));%P(xd) for each bin
    end
    P(j+1,i)=sum(intj);% P(xd)
    
end
%Ensure there are no small negative values
for j=1:no_b+1
    pj=P(j,:);
    pj(find(pj<0))=0;
    P(j,:)=pj;
end

nc=trapz(xds,P(end,:));%normalisation constant
P=P/nc;%normalise