size = 180;
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
    if (A(i,1) > 0) && (A(i,1) <= 0.15)
        A(i,1) = 0;
    end
    if (A(i,1) < 0) && (A(i,1) >= -0.15)
        A(i,1) = 0;
    end
    if (A(i,2) > 0) && (A(i,2) <= 0.15)
        A(i,2) = 0;
    end
    if (A(i,2) < 0) && (A(i,2) >= -0.15)
        A(i,2) = 0;
    end
    if (A(i,3) > 0) && (A(i,3) <= 0.15)
        A(i,3) = 0;
    end
    if (A(i,3) < 0) && (A(i,3) >= -0.15)
        A(i,3) = 0;
    end
    
    Awc = p1*A(i,:)';
    u(i,:) = u(i-1,:)+Awc'*dt;
    pos(i,:) = pos(i-1,:);
    pos(i,:) = pos(i,:) + Awc'*(dt*dt)/2 + u(i-1,:)*dt ;
    posn(i) = norm(pos(i,:));
    An(i) = norm(A(i,:));
    un(i) = norm(u(i,:));
    Gguess(i,:) = p1(1,:);
end
plot(1:size,posn);
% M = [M Gguess];
% yy = fft(A);