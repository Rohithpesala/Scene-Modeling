a = dlmread('LACC_20160418_150028_786.txt',',');
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
for i=2:length(a)
    dt = arr(i)/10^9;
    v(i,:) = v(i-1,:)+tau*a(i,1:3)*dt;
    x(i,:) = x(i-1,:)+v(i-1,:)*dt+0.5*a(i,1:3)*(dt)^2;
end
xn = zeros(length(a),3);
for i=1:length(x)
    xn(i) = norm(x(i));
end
plot(xn);
% sum1 = 0;
% for i = sp:ep
%     ix = x(round((i+interval)*scale),:)-x(round(i*scale),:);    %row vector
%     vx = Ot(:,i);                   %column vector
%     sum1 = sum1 + vx'*ix';
% end
% lam = sum1/(ep-sp+1);