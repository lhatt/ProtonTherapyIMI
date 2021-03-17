function a=sum_stat(k)
% Returns the summary statistic with length k- it is not normalised
a = [1; 2*rand(k,1)-1];%Range [-1,1]; 
a = a'.*(linspace(1,k+2,k+1).^(-1));
a=a(2:end);