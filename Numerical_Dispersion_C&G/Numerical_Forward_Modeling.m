% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% 2D acoustic wave propagation                                             %
% Layered Velocity                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%  Modified from Xin Wang (2012)    -- Center for Subsurface Imaging and Fluid Modeling (CSIM),
%                                      King Abdullah University of Science and Technology, All rights reserved.

%  purpose:  2DTDFD solution to acoustic wave equation with accuracy of 2-8
%            use the absorbing boundary condition
%
%  IN   v(:,:) -- velocity,      nbc         -- grid number of boundary
%       dx     -- grid intervel, nt          -- number of sample
%       dt     -- time interval, s(:)        -- wavelet
%       sx,sz  -- src position,  gx(:),gz(:) -- rec position
%       isFS   -- if is Free Surface condition

clc
clear all

% global vel;
% global x;   global z;
% global t;
% global ng;  global g;
% global nx;  global nz;

% Kecepatan Berlapis
vel=[repmat(2200,[1,55]), repmat(2700,[1,55]), repmat(3500,[1,55]), repmat(4100,[1,55]), repmat(4700,[1,55]), repmat(5300,[1,46])];
% vel=[repmat(4200,[1,55]), repmat(4700,[1,55]), repmat(5300,[1,55]), repmat(5800,[1,55]), repmat(6300,[1,55]), repmat(6900,[1,46])];
vel=repmat(vel',[1 451]);

% vel=[repmat(2200,[1,55]), repmat(2400,[1,55]), repmat(2600,[1,55]), repmat(2900,[1,55]), repmat(3200,[1,55]), repmat(3400,[1,46])];
% vel=repmat(vel',[1 451]);

% % Kecepatan Berlapis (3-Lapis)
% vel=[repmat(2700,[1,55]), repmat(2700,[1,55]), repmat(4100,[1,55]), repmat(4100,[1,55]), repmat(5300,[1,55]), repmat(5300,[1,46])];
% % vel=[repmat(2200,[1,55]), repmat(2200,[1,55]), repmat(2700,[1,55]), repmat(2700,[1,55]), repmat(3500,[1,55]), repmat(3500,[1,46])];
% vel=repmat(vel',[1 451]);

% % Kecepatan Homogen
% vel=[repmat(2700,[1,55]), repmat(2700,[1,55]), repmat(2700,[1,55]), repmat(2700,[1,55]), repmat(2700,[1,55]), repmat(2700,[1,46])];
% vel=repmat(vel',[1 451]);

[nz,nx]=size(vel); 
nbc = 80;
% load velocityModel
% [nz,nx] = size(velocityModel);

M = 5;
% dx = 7.625; 
dx = 15.625; 
% dx = 20.625; 
x = (0:nx-1)*dx; 

% dz = 7.625;   
dz = 15.625;
% dz = 20.625;
z = (0:nz-1)*dz;

sx = ((nx-1)/2)*dx; 
% sz = 650;
sz = 1300;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

isx=round(sx/dx)+1+nbc;
isz=round(sz/dx)+1+nbc;

% gx = 0:dx:(nx-1)*dx; 
gx = 0:dx:(nx-1)*dx;
gz = zeros(size(gx)); 
ng = numel(gx); 
g = 1:ng; 

figure(1);
set(gcf,'position',[0 0 1000 400]);
imagesc(x,z,vel);
set(gca,'FontSize',16,'fontWeight','bold')
% imagesc(x,z,velocityModel);
colorbar;
xlabel('X (m)'); ylabel('Z (m)'); 
% title('Homogeneous Velocity Model \Deltax = 15.625 m Sz = 1300 m');
title({'Layered Velocity Model (6-Layers) \Deltax = 15.625 m Sz = 1300 m'});

hold on
plot(gx,gz,'r^',sx,sz,'rp');

% Ricker Source
nt = 15001;
% nt = 30001;
dt = 0.0002;
% dt = 0.0001;
t = (0:nt-1)*dt;
var = vel*dt/dx;
isFS = false;
freq = 25; 
% freq = 30; 
s=ricker(freq,dt,nt);
% s=ricker2(freq,nt,dt,t,0);

if var <= sqrt(3/8)
    disp('Lanjutkan                           [  OK  ]')
elseif var > sqrt(3/8) 
    disp('Berhenti                           [ GAGAL ]'), return
end

seis=zeros(nt,numel(gx));
% [nz,nx]=size(v);

% c1 = -2; c2 = 1;
% c1 = -205.0/72.0; c2 = 8.0/5.0; c3 = -1.0/5.0; c4 = 8.0/315.0; c5 = -1.0/560.0;

% global M;
% M = 5;
[c,c1]=linear_eq(M);
% [cd1,c1]=linear_eq3(dt,dx,v);
% [cd1,c1]=linear_eq31;
% [cd1,c1]=linear_eq32;
% [cd1,c1]=linear_eq4;

% c2 = c(1);
% c2 = cd1(:,:,1);

% c3 = c(2);
% c3 = cd1(:,:,2);

% c4 = c(3); 
% c4 = cd1(:,:,3);

% c5 = c(4);
% % c5 = cd1(:,:,4);

% c6 = c(5);
% % c6 = cd1(:,:,5);

% c7 = c(6);
% % c7 = cd1(:,:,6);

% c8 = c(7);
% % c8 = cd1(:,:,7);

% c9 = c(8);
% % c9 = cd1(:,:,8);

% c10 = c(9);
% % c10 = cd1(:,:,9);

% c11 = c(10);
% % c11 = cd1(:,:,10);

% c12 = c(11);
% % c12 = cd1(:,:,11);

% c13 = c(12);
% % c13 = cd1(:,:,12);

% c14 = c(13);
% % c14 = cd1(:,:,13);

% setup ABC and temperary variables
vel = padvel(vel,nbc); 
% c1=padvel(c1,nbc);  c2=padvel(c2,nbc);      c3=padvel(c3,nbc);      c4=padvel(c4,nbc);
% c5=padvel(c5,nbc);  c6=padvel(c6,nbc);      c7=padvel(c7,nbc);      c8=padvel(c8,nbc);
% c9=padvel(c9,nbc);  c10=padvel(c10,nbc);    c11=padvel(c11,nbc);    c12=padvel(c12,nbc);
% c13=padvel(c13,nbc);  c14=padvel(c14,nbc);
% c1=padvel(c1,nbc);      cd1 = padvel(cd1(:,:,1:M),nbc);

abc=AbcCoef2D(vel,nbc,dx);
alpha=(vel*dt/dx).^2; 
kappa=abc*dt;
temp1=2+2*c1.*alpha-kappa; 
temp2=1-kappa;
beta_dt = (vel*dt).^2;
s=expand_source(s,nt);
[isx,isz,igx,igz]=adjust_sr(sx,sz,gx,gz,dx,nbc);
p1=zeros(size(vel)); p0=zeros(size(vel));

% Time Looping
% vidObj1 = VideoWriter('Forward Modeling of Homogeneous Model dx = 15.625 m Freq = 25 Hz 4th Order Sz = 1300 m.avi');
% open(vidObj1)
vidObj1 = VideoWriter('Forward Modelling of Layered Model (6-Layers) dx = 15.625 m Freq = 25 Hz 10th Order Sz = 1300 m.avi');
open(vidObj1);
for it=1:nt
    sum = 0;
    for m = 1:M
        sum = sum + c(m).*(circshift(p1,[0,m,0])+circshift(p1,[0,-m,0])+circshift(p1,[m,0,0])+circshift(p1,[-m,0,0]));
        % sum = sum + cd1(:,:,m).*(circshift(p1,[0,m,0])+circshift(p1,[0,-m,0])+circshift(p1,[m,0,0])+circshift(p1,[-m,0,0]));
        p = temp1.*p1 - temp2.*p0 + alpha.*(sum);
    end

    p(isz,isx) = p(isz,isx) + beta_dt(isz,isx) * s(it);
    % dipole source
    
    %p(isz-2,isx) = p(isz-2,isx) - beta_dt(isz-2,isx) *wavelet(it);
    if isFS
        p(nbc+1,:)=0.0;
        p(nbc:-1:nbc-3,:) = - p(nbc+2:nbc+5,:);
    end
    for ig=1:ng
        seis(it,ig)=p(igz(ig),igx(ig));
    end
    if mod(it,100)==1
        figure(2);
        % colormap(gray);
        set(gcf,'position',[0 0 1200 1000]);
        imagesc(x,z,p(nbc+1:nbc+nz,nbc+1:nbc+nx));
        set(gca,'FontSize',14,'fontWeight','bold');
        % figure_title=['Wave Propagation of Homogeneous Model \Deltax = 15.625 m Freq = 25 Hz 4th Order, Sz = 1300 m, t=', num2str((it-1)*dt),'s'];
        figure_title = ['Wave Propagation of Layered Model (6-layers) \Deltax = 15.625 m Freq = 25 Hz, 10th Order, Sz = 1300 m, t=',num2str((it-1)*dt),' s'];
        title(figure_title);
        % title({'Wave Propagation of Layered Model (6-layer)';' \Deltax = 15.625 m Freq = 25 Hz Orde-10 Sz = 1300 m, t=',num2str((it-1)*dt),' s'});
        ylabel('Z (m)');
        xlabel('X (m)');
        caxis([-0.25 0.25]);
        
        hold on
        colorbar;
        % set(gca,'CLim',[-10^-9 10^-9])
        % plot(gx,gz,'r^',x(isx),z(isz),'rp');
        plot(gx,gz,'r^',sx,sz,'rp');
        writeVideo(vidObj1,getframe(gcf));
        
        figure(3);
        imagesc(g,t,seis);
        set(gca,'FontSize',14,'fontWeight','bold');
        % figure_title=['Seismic Profie of Homogeneous Model \Deltax = 15.625 m Freq = 25 Hz 4th Order, Sz = 1300 m, t=', num2str((it-1)*dt),' s'];
        figure_title=['Seismic Profie of Layered Model (6-layers) \Deltax = 15.625 m Freq = 25 Hz, 10th Order, Sz = 1300 m, t=',num2str((it-1)*dt),' s'];
        title(figure_title);
        % title({'Seismic Profie of Layered Model (6-layer)';' \Deltax = 15.625 m Freq = 25 Hz Orde-10 Sz = 1300 m t=',num2str((it-1)*dt),' s'});
        ylabel('Time (s)');
        xlabel('Offset (m)');
        caxis([-0.25 0.25]);
        pause(0.2);
    end
    p0=p1;
    p1=p;
end
close(vidObj1)

figure(4);
diseis(seis,2,g,t)
set(gca,'FontSize',15,'fontWeight','bold');
% title({'Seismic Profile 0f Homogeneous Model - Wiggle';' \Deltax = 15.625 m Freq = 25 Hz 4th Order Sz = 1300 m'});
title({'Seismic Profile of Layered Model (6-Layers) - Wiggle';' \Deltax = 15.625 m Freq = 25 Hz, 10th Order Sz = 1300 m'});
xlabel('g #');
ylabel('Time (s)');

ngx = floor(ng/2);
figure(5);
plot(t,seis(:,ngx));
set(gca,'FontSize',15,'fontWeight','bold');
% title({'Single Seismogram 0f Homogeneous Model ';' \Deltax = 15.625 m Freq = 25 Hz 4th Order Sz = 1300 m'});
title({'Single Seismogram of Layered Model (6-Layers)';' \Deltax = 15.625 m Freq = 25 Hz, 10th Order Sz = 1300 m'});
xlabel('t (s)');
ylabel('Amplitude');








