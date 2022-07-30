close all
clear
clc

%% Parameters
h=0.01;             % Spatial Step (dx=dy=h)
L=1;                % Length
n=L/h+1;            % Number of Steps
x=0:h:L;            % X Bound
y=0:h:L;            % Y Bound
U=1;                % Lid Velocity
dt=0.002;           % Temporal Step
endTime=3;          % Final Time
Re=500;             % Reynolds Number
nu=(U*L)/Re;        % Courant Number
% Initialize 
omega=zeros(n,n);   % Initialize Omega
psi=zeros(n,n);     % Initialize Stream Function
u=zeros(n,n);       % X Direction Velocity
v=zeros(n,n);       % Y Direction Velocity
s=0;                % Time Loop's Counter
t=0;                % Time
iter=0;             % Iteration
q=0.5;              % blending factor
%% initial Condition
for i=1:n
    for j=1:n
        psi(i,j,1)=0;
        omega(i,j,1)=0;
        if j==n
            psi(i,j,1)=0;
            omega(i,j,1)=-2*U/h;
        end
    end
end
%% time loop
while t<endTime
    vortold=omega;
    psiold=psi;
    iter=iter+1;
    t=t+dt;
    s=s+1;
    %vorticity
    for i=1:n
        for j=1:n
            if j~=1 && j~=n && i~=1 && i~=n
central_ypsi=(psi(i,j+1)-psi(i,j-1));
umin=min(central_ypsi,0);
uplus=max(central_ypsi,0);
central_xpsi=(psi(i+1,j)-psi(i-1,j));
vmin=min(central_xpsi,0);
vplus=max(central_xpsi,0);
                                                       % 1st-O upwind for blending
back_xomega= (omega(i,j)-omega(i-1,j))/h;
forward_xomega = (omega(i+1,j)-omega(i,j))/h;
back_yomega= (omega(i,j)-omega(i,j-1))/h;
forward_yomega = (omega(i,j+1)-omega(i,j))/h;
omega(i,j)=omega(i,j)+dt*(-central_ypsi*(omega(i+1,j)-omega(i-1,j))/(4*h^2)+q*(uplus*back_xomega + umin*forward_xomega)+central_xpsi*(omega(i,j+1)-omega(i,j-1))/(4*h^2)+q*(vplus*back_yomega + vmin*forward_yomega)+nu*(omega(i+1,j)+omega(i,j+1)+omega(i-1,j)+omega(i,j-1)-4*omega(i,j))/h^2);
            elseif j==1 && i~=1 && i~=n
                omega(i,j)=2*(psi(i,j)-psi(i,j+1))/h^2;          % Bottom Condition
            elseif j==n && i~=1 && i~=n
                omega(i,j)=2*(psi(i,j)-psi(i,j-1))/h^2-2*U/h;    % Top Condition
            elseif i==1 && j~=1 && j~=n
                omega(i,j)=2*(psi(i,j)-psi(i+1,j))/h^2;          % Left Condition
            elseif i==n && j~=1 && j~=n
                omega(i,j)=2*(psi(i,j)-psi(i-1,j))/h^2;          % Right Condition
            elseif (i==1 && j==1)                                % Set corners equal to neighbors or zero
                omega(i,j)=0;
            elseif (i==n && j==1)
                omega(i,j)=0;
            elseif (i==1 && j==n)
                omega(i,j)=(omega(i+1,j));
            elseif (i==n && j==n)
                omega(i,j)=(omega(i-1,j));
            end
        end
    end
    %streamfunction 
    for i=1:n
        for j=1:n
            if j~=1 && j~=n && i~=1 && i~=n
                psi(i,j)=0.25*(h^2*omega(i,j)+psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1));
            else
                psi(i,j)=0;
            end
        end
    end
    %calculate velocities
    for i=1:n
        for j=1:n
            if j~=1 && j~=n && i~=1 && i~=n
                u(i,j)=(psi(i,j+1)-psi(i,j-1))/(2*h);
                v(i,j)=(psi(i+1,j)-psi(i-1,j))/(2*h);
            elseif j==n
                u(i,j)=1;
                v(i,j)=0;
            else
                u(i,j)=0;
                v(i,j)=0;
            end
        end
    end
vortout{s}(:,:)=omega;
uout{s}(:,:)=u;
vout{s}(:,:)=v;
end
%% plots
subplot(2,2,1)
quiver(x,y,u,v)
xlabel('x');
ylabel('y');
title(sprintf('Velocity Field Re=500 q=%0.2f t=%0.2f',q,t));

subplot(2,2,2)
contourf(x,y,psi',20,'LineColor','none' )
xlabel('x');
ylabel('y');
colorbar
title(sprintf('Streamlines Re=500 q=%0.2f t=%0.2f',q,t));

subplot(2,2,3)
contourf(x,y,abs(omega'),40,'LineColor','none' )
xlabel('x');ylabel('y');
colorbar
title(sprintf('Vorticity Re=500 q=%0.2f t=%0.2f',q,t));