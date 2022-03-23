function [c,c1]=linear_eq(M)

% Author : Syamsurizal Rizal

% Purpose     : Determining finite difference coefficient of 2D TDFD
%               acoustic wave equation


% M = 5;
% i = 2:M;
for i = 2:M
    for j = 1:M
        A(1,j)=j^2;
        A(i,j)=A(1,j).^i;
        % A(i+1,j)=A(1,j);
    end  
    % A(i+1,j)=A(1,j);
    R(i)=0;
    % R(i)=r.^i;
end
% A
R(1)=1;
% R'
c = inv(A)*R';
fprintf('c = %7.30f \n',c);
c1=-2*sum(c);
fprintf('c1 = %7.10f \n',c1);
