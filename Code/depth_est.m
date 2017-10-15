fno = 13;
lam = tlam(3);
P = OP(:,:,fno);
P(:,4) = P(:,4)*lam;
m3d = [];
for i=1:length(l_p)
    t_p = ch_p(i,:)';
    t_p1 = l_p(i,:)';
    temp1 = pinv(K*P)*t_p;
    temp0 = pinv(K*[eye(3) [0;0;0]])*t_p1;
    A = [temp0(1:3) -temp1(1:3)];
    B = P(1:3,4);
    ot1 = pinv(A)*B;
    o1a = ot1(1)*temp0(1:4);
    if inl(i) == 1
        m3d = [m3d; o1a(1:3)'];
    end
end
m3d = [m3d ones(length(m3d),1)];
% figure
% plot(m3d(:,1),m3d(:,3),'.')
m2d = zeros(length(m3d),2);
for i=1:length(m3d)
    temp = (K*P)*m3d(i,:)';
    temp = temp/temp(3);
    m2d(i,:) = round(temp(1:2));
end
