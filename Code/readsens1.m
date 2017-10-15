a = dlmread('LACC_20160528_101022_492.txt',',');
arr = zeros(length(a),1);
tau = 0.2;
scale = length(a)/nf;
for i=2:length(a)
    arr(i) = a(i,4)-a(i-1,4);
%     arr(i) = norm(a(i,1:3));
end
iter = 1:length(a);
% plot(iter,a(:,3));
v = zeros(length(a),3);
x = zeros(length(a),3);
%filter
aa = 1;
bb = ones(20,1)*0.05;
f = filter(bb,aa,a(:,1:3));
[xx,yy] = butter(20,0.2,'low');
f = filter(xx,yy,a(:,1:3));
for i = 1:length(f)
    ff(i) = norm(f(i,:));
end
a = f;
%basic model
for i=2:length(a)
    dt = arr(i)/10^9;
    v(i,:) = v(i-1,:)+a(i-1,1:3)*dt;
    x(i,:) = x(i-1,:)+v(i-1,:)*dt+0.5*a(i-1,1:3)*(dt)^2;
end
%verlet model
% for i=3:length(a)
%     dt = arr(i)/10^9;
%     x(i,:) = 2*x(i-1,:)-x(i-2,:)+a(i-1,1:3)*(dt)^2;
% end

xn = zeros(length(a),3);
an = zeros(length(a),3);
for i=1:length(x)
    xn(i) = norm(x(i));
    an(i) = norm(a(i,1:3));
end
% plot(xn);
sum1 = 0;
nvx = [];
nix = [];
sp = 28;
ep = 30;
for i = sp:ep
    ix = x(round((i+interval)*scale),:)-x(round(i*scale),:);    %row vector
    vx = Ot(:,i);                   %column vector
    sum1 = sum1 + vx'*ix';
    nvx = [nvx; vx'];
    nix = [nix; ix];
end
lam = sum1/(ep-sp+1);
% figure
% plot(an); hold on
nvx = lam*nvx;
v1 = var(nvx);
v2 = var(nix);
k1 = v1(1)*v2(1)/(v1(1)+v2(1));
k2 = v1(2)*v2(2)/(v1(2)+v2(2));
k3 = v1(3)*v2(3)/(v1(3)+v2(3));
xf = zeros(size(nvx));
for i = 1:(ep-sp+1)
    xf(i,:) = nvx(i,:)./v1 + nix(i,:)./v2;
    xf(i,:) = xf(i,:).*[k1 k2 k3];
end
txf = txf + sum(xf)
tlam = [tlam lam]
% plot(x(:,1),'r')