clc;
clear all;

% Author      : Syamsurizal Rizal

% vel         : Velocity
% dt          : time step
% freq        : frequency


vel = 3000;
dt = 0.00015;
% dt = 0.001;
dx = 1.0:0.25:50;
r = vel*dt./dx;

freq = 25;
omega = 2*pi*freq;
k1 = omega./vel;
r1 = r';
B = (cos(omega*dt)-1)./(r1.^2);
theta = pi/4;

%% M1

M1 = 2;
% M1 = 9;
ca1=zeros(M1,length(dx));
cb1=zeros(length(dx),M1);
for k = 1:length(dx)
    for i = 2:M1
        for j = 1:M1
            A1 = cos(theta)^2;    A2 = sin(theta)^2;
            A(1,j) = (j^2)*(A1+A2);
            A(i,j) = (A(1,j).^i).*(A1.^i + A2.^i);
        end
        R(i)=r(k).^(2*i-2);
        % R(i)=0;
    end  
        
    R(1)=1;
    c = inv(A)*R'
    fprintf('c = %7.30f \n',c);
    
    ca1(:,k) = c(:,1);
    cb1(k,:) = ca1(:,k)';
       
    c11(k)=-2*sum(c(:,1));
end
d1(:,:,1)=cb1;

S11 = cos(theta);   S12 = sin(theta);
m1 = 1:M1;
S1 = -2+cos(m1.*k1.*dx'.*S11)+cos(m1.*k1.*dx'.*S12);

% sum =0;
sum1=zeros(length(dx),1);
for m1 = 1:M1
    sum1 = sum1 + d1(:,m1).*S1(:,m1);
end
error1 = (sum1-B)./(dx'.^2);

%% M2

M2 = 4;
% M2 = 10;
ca2=zeros(M2,length(dx));
cb2=zeros(length(dx),M2);
for k = 1:length(dx)
    for i = 2:M2
        for j = 1:M2
            A1 = cos(theta)^2;    A2 = sin(theta)^2;
            A(1,j) = (j^2)*(A1+A2);
            A(i,j) = (A(1,j).^i).*(A1.^i + A2.^i);
        end
        R(i)=r(k).^(2*i-2);
        % R(i)=0;
    end  
        
    R(1)=1;
    c = inv(A)*R'
    fprintf('c = %7.30f \n',c);
    
    ca2(:,k) = c(:,1);
    cb2(k,:) = ca2(:,k)';
       
    c12(k)=-2*sum(c(:,1));
end
d2(:,:,1)=cb2;

m2 = 1:M2;
S2 = -2+cos(m2.*k1.*dx'.*S11)+cos(m2.*k1.*dx'.*S12);

% sum =0;
sum2=zeros(length(dx),1);
for m2 = 1:M2
    sum2 = sum2 + d2(:,m2).*S2(:,m2);
end
error2 = (sum2-B)./(dx'.^2);

%% M3

M3 = 6;
% M3 = 11;
ca3=zeros(M3,length(dx));
cb3=zeros(length(dx),M3);
for k = 1:length(dx)
    for i = 2:M3
        for j = 1:M3
            A1 = cos(theta)^2;    A2 = sin(theta)^2;
            A(1,j) = (j^2)*(A1+A2);
            A(i,j) = (A(1,j).^i).*(A1.^i + A2.^i);
        end
        R(i)=r(k).^(2*i-2);
        % R(i)=0;
    end  
        
    R(1)=1;
    c = inv(A)*R'
    fprintf('c = %7.30f \n',c);
    
    ca3(:,k) = c(:,1);
    cb3(k,:) = ca3(:,k)';
       
    c13(k)=-2*sum(c(:,1));
end
d3(:,:,1)=cb3;
  
m3 = 1:M3;
S3 = -2+cos(m3.*k1.*dx'.*S11)+cos(m3.*k1.*dx'.*S12);

% sum =0;
sum3=zeros(length(dx),1);
for m3 = 1:M3
    sum3 = sum3 + d3(:,m3).*S3(:,m3);
end
error3 = (sum3-B)./(dx'.^2);

%% M4

M4 = 8;
% M4 = 12;
ca4=zeros(M4,length(dx));
cb4=zeros(length(dx),M4);
for k = 1:length(dx)
    for i = 2:M4
        for j = 1:M4
            A1 = cos(theta)^2;    A2 = sin(theta)^2;
            A(1,j) = (j^2)*(A1+A2);
            A(i,j) = (A(1,j).^i).*(A1.^i + A2.^i);
        end
        R(i)=r(k).^(2*i-2);
        % R(i)=0;
    end  
        
    R(1)=1;
    c = inv(A)*R'
    fprintf('c = %7.30f \n',c);
    
    ca4(:,k) = c(:,1);
    cb4(k,:) = ca4(:,k)';
       
    c14(k)=-2*sum(c(:,1));
end
d4(:,:,1)=cb4;

m4 = 1:M4;
S4 = -2+cos(m4.*k1.*dx'.*S11)+cos(m4.*k1.*dx'.*S12);

% sum =0;
sum4=zeros(length(dx),1);
for m4 = 1:M4
    sum4 = sum4 + d4(:,m4).*S4(:,m4);
    % B = (cos(omega*dt)-1)./r1^2;
end
error4 = (sum4-B)./(dx'.^2);

%% Plot Grafik

plot(dx,error1,dx,error2,dx,error3,dx,error4,'linewidth',2);
set(gca,'FontSize',14,'fontWeight','bold')
xlabel('\Deltah(m)');
ylabel('Numerical Error');
title({'Relation Between Numerical Errors';'To Grid Size';'Velocity = 3000 m/s, \Deltat = 0,00015s, Frequency = 25 Hz, \Theta = \pi/4'});
legend({'M = 2','M = 4','M = 6','M = 8'},'Location','northwest');
% legend({'M = 9','M = 10','M = 11','M = 12'},'Location','northwest');

% legend({'M=9','M=10','M=11','M=12'},'Location','northwest','NumColumns','2');
% legend({'M = 9','M = 10','M = 11','M =1 2'},'Location','northwest');
