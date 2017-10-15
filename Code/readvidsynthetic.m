newnp = 12;            %%% number of points to use for correspondences
%% Calibrated intrinsic params
K = [2088 0 1269; 0 2066 969;0 0 1];

loc1 = pp1';
cor_loc = pp2';

A = zeros(newnp,9);
for j = 1:newnp
    ai = cor_loc(j,:)'*loc1(j,:);
    %ai = ai';
    ai = ai(:);
    ai = ai';
    A(j,:) = ai;
end

%% Finding the fundamental matrix as the eigen vector corresponding to the least eigen value
%     [U,S,V] = svd(A);
[Ve,De] = eig(A'*A);
f = Ve(:,1)';
%     f = V(:,9)';

F = zeros(3,3);
F(1,:) = f(1:3);
F(2,:) = f(4:6);
F(3,:) = f(7:9);
%     F = T_1'*F*T_2;
[U,S,V] = svd(F');
S(3,3) = 0;         % equating the last singular value to 0 as F is rank 2
Ft = U*S*V';

% [fRANSAC, inl] = estimateFundamentalMatrix(loc1(:,1:2),cor_loc(:,1:2));
% Ft(:,:,i-1) = fRANSAC;

%%%Tranforming back the F matrix to scale up to remove normalization
% Ft(:,:,i-1) = T_1*Ft(:,:,i-1)*T_2';
E = K'*Ft*K;   % Converting to essential matrix

[U,S,V] = svd(E);
S(1,1) = 1;
S(2,2) = 1;
S(3,3) = 0;
E = U*S*V';

%% Getting the rotation and translation vectors from essential matrix
W = [0 -1 0;1 0 0;0 0 1];
t1 = U(:,3);
t2 = -U(:,3);
V(:,3) = -V(:,3);
R1 = U*W*V';
R2 = U*W'*V';
P1 = [R1 t1];
P2 = [R2 t1];
P3 = [R1 t2];
P4 = [R2 t2];
sc = zeros(4,newnp);

%% validating the obtained projection matrices
loc1(1,:) = [440  950    1.0000]; %experiment
for i = 1:1%newnp
    t_p = cor_loc(i,:)';
    t_p1 = loc1(i,:)';
    temp1 = pinv(K*P1)*t_p;
    temp0 = pinv(K*[eye(3) [0;0;0]])*t_p1;
    A = [temp0(1:3) -temp1(1:3)];
    B = P1(1:3,4);
    ot1 = pinv(A)*B;
    o1a = ot1(1)*temp0(1:4);
    o1b = K*P1*o1a;
%     sc(1,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
    temp2 = pinv(K*P2)*t_p;
    A = [temp0(1:3) -temp2(1:3)];
    B = P2(1:3,4);
    ot2 = pinv(A)*B;
    o2a = ot2(1)*temp0(1:4);
    o2b = K*P2*o2a;
%     sc(2,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
    temp3 = pinv(K*P3)*t_p;
    A = [temp0(1:3) -temp3(1:3)];
    B = P3(1:3,4);
    ot3 = pinv(A)*B;
    o3a = ot3(1)*temp0(1:4);
    o3b = K*P3*o3a;
%     sc(3,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
    temp4 = pinv(K*P4)*t_p;
    A = [temp0(1:3) -temp4(1:3)];
    B = P4(1:3,4);
    ot4 = pinv(A)*B;
    o4a = ot4(1)*temp0(1:4);
    o4b = K*P4*o4a;
%     sc(4,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
end

