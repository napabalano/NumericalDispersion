function [cd1,c1]=linear_eq3(dt,dx,vel)

% Author      : Syamsurizal Rizal

% Purpose     : Determining finite diffference coefficient of 2D TDFD
%               acoustic wave equation

% global M;

% vel=[repmat(2200,[1,55]), repmat(2700,[1,55]), repmat(3500,[1,55]), repmat(4100,[1,55]), repmat(4700,[1,55]), repmat(5300,[1,46])];
% vel=repmat(vel',[1 451]);
[N1,N2] = size(vel);
B = zeros(N1*N2,1);
for i=1:N1
    for j=1:N2
        B((i-1)*N2 + j) = vel(i,j);
    end
end
% dt = 0.0002;  dx = 15.625;
r = B*dt/dx;

M = 12; 
c = zeros(M,1);
for k = 1:length(B)
    for i = 2:M
        for j = 1:M
            A(1,j) = (j)^2;
            A(i,j)=A(1,j).^i;
        end
        R(i)=r(k).^(2*i-2);
    end  
        
    R(1)=1;
    c=inv(A)*R'
    fprintf('c = %7.30f \n',c);
    
    ca(:,k) = c(:,1);
    cb(k,:) = ca(:,k)';
    
    c1(k)=-2*sum(c(:,1));
end
c1 = reshape(c1(:),N1,N2);
for i = 1:M
    cd1(:,:,i) = reshape(cb(:,i),N1,N2);
end
