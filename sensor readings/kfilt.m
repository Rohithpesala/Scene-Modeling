size = 81;
M = csvread('idle.csv',1,3,[1,3,size,5]);
M(:,3) = -1*M(:,3);
G = csvread('idle.csv',1,9,[1,9,size,11]);
G(:,3) = -1*G(:,3);
A = csvread('idle.csv',1,6,[1,6,size,8]);
A(:,3) = -1*A(:,3);
dt = 0.1;
V = M(1,:);
V = V/norm(V);
rndv = [1,1,1];
Gguess = zeros(size,3);
h1 = cross(rndv,V);
h1 = h1/norm(h1);
h2 = cross(V,h1);
h2 = h2/norm(h2);
p1 = [V;h1;h2];

u = zeros(size,3);
pos = zeros(size,3);
posn = zeros(size,1);
An = zeros(size,1);
un = zeros(size,1);
F = zeros(9,9);
F(1,:) = [1 0 0 dt 0 0 dt^2/2 0 0];
F(2,:) = [0 1 0 0 dt 0 0 dt^2/2 0];
F(3,:) = [0 0 1 0 0 dt 0 0 dt^2/2];
F(4,:) = [0 0 0 1 0 0 dt 0 0];
F(5,:) = [0 0 0 0 1 0 0 dt 0];
F(6,:) = [0 0 0 0 0 1 0 0 dt];
F(7,:) = [0 0 0 0 0 0 1 0 0];
F(8,:) = [0 0 0 0 0 0 0 1 0];
F(9,:) = [0 0 0 0 0 0 0 0 1];
X_old = zeros(9,1);
R = [0.2 0 0;0 0.2 0;0 0 0.2];
H = [0 0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 1];
P_old = diag(zeros(9,1)+0.2);
pos(1,:) = [0 0 0];

for i=2:size
    scale_n = 10;
    phi_x = (G(i,1)*dt)/180*(pi/scale_n);
    phi_y = (G(i,2)*dt)/180*(pi/scale_n);
    phi_z = (G(i,3)*dt)/180*(pi/scale_n);
    R_x = [1 0 0;0 cos(phi_x) -sin(phi_x);0 sin(phi_x) cos(phi_x)];
    R_y = [cos(phi_y) 0 sin(phi_y);0 1 0;-sin(phi_y) 0 cos(phi_y)];
    R_z = [cos(phi_z) -sin(phi_z) 0;sin(phi_z) cos(phi_z) 0;0 0 1];
    Rt = R_x*R_y*R_z;
    Rf = Rt^scale_n;
    p1 = p1*inv(Rf);
    
%     Awc = p1*A(i,:)';
    Awc = A(i,:)';
    %% Predict equations
    X_pred = F*X_old;
    P_pred = F*P_old*F';
    %% Kalman gain
    Kt = P_pred*H'*inv(H*P_pred*H' + R);
    %% Update equations
    X_upd = X_pred + Kt*(Awc - H*X_pred);
    P_upd = P_pred - Kt*H*P_pred;
    %% reassigning vars
    X_old = X_upd;
    P_old = P_upd;
        
%     u(i,:) = u(i-1,:)+Awc'*dt;
    pos(i,:) = X_old(1:3);
    posn(i) = norm(pos(i,:));
%     pos(i,:) = pos(i,:) + Awc'*(dt*dt)/2 + u(i-1,:)*dt ;
%     posn(i) = norm(pos(i,:));
%     An(i) = norm(A(i,:));
%     un(i) = norm(u(i,:));
    Gguess(i,:) = p1(1,:);
end
plot(1:size,posn);
% M = [M Gguess];
% yy = fft(A);