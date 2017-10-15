v = VideoReader('VID_20160407_140014_798.mp4');
nf = v.NumberOfFrames;
h = v.Height;
w = v.Width;
np = 300;            %%% number of points to use for correspondences
Ft = zeros(3,3,nf-1);
tfactor = 2;      %%% tolerance factor
rem = 0;
sfactor = 6;
cfactor = 5;

%% Calibrated intrinsic params
% K = [2088 0 1269; 0 2066 969;0 0 1];
K = [4212 0 2564; 0 4169 1573;0 0 1];

tic
%%% looping over each frame
lp = 64;
for i = lp:lp
    i
    
    %%% two frames from video
    v1 = im2double(rgb2gray(read(v,i-1)));
    v2 = im2double(rgb2gray(read(v,i+20)));
    vv2 = v2;
%     tv1 = rgb2gray(read(v,i));
%     vz = zeros(1,w);
%     tv2 = [vz ; tv1(1:h-1,:)];
%     imwrite(v1,'initframe.jpg','jpg');
%     imwrite(v2,'finframe.jpg','jpg');
%     v1 = imread('initframe.jpg','jpg');
%     v2 = imread('finframe.jpg','jpg');
%     imshow(v1);
%     figure
%     imshow(v2);
    
    %% Finding feature points in first image
    points1 = detectSURFFeatures(v1);
    points22 = detectSURFFeatures(vv2);
    p1 = points1.selectStrongest(np);
    loc1 = zeros(np,2);
    met = zeros(np,1);
    scrs = zeros(np*3,np);
    cor_loc = zeros(np,2);
    for j = 1:np
        loc1(j,:) = round(p1(j).Location);
        met(j) = round(p1(j).Scale);
    end
    met = met*sfactor;
    
    %% adding padding for second frame with max radius to avoid conflicts. Please neglect this padding
    met_m = ceil(max(met))*2;
    vz = zeros(met_m,w);
    v2 = [vz ; v2 ; vz];
    hz = zeros(h+2*met_m,met_m);
    v2 = [hz v2 hz];
    
    %% Findiing feature points in second frame
    points2 = detectSURFFeatures(v2);
    p2 = points2.selectStrongest(np*3);
    loc2 = zeros(np*3,2);     
    for j = 1:np*3
        loc2(j,:) = round(p2(j).Location);        
    end
    
    %% Finding the correspondences by comparing two frames feature points
    ttemp = 0
    rem1 = 0
    for jj = 1:np
        j = jj-rem1;
        if (loc1(j,2)-met(j) >0  && loc1(j,1)-met(j)>0 && loc1(j,2)+met(j) <=h && loc1(j,1)+met(j)<=w)
            sec1 = v1(loc1(j,2)-met(j) : loc1(j,2)+met(j), loc1(j,1)-met(j) : loc1(j,1)+met(j));
            ms = inf;
            mind = 1;
            for k = 1:np*3
                y1 = loc2(k,2)-met(j);
                y2 = loc2(k,2)+met(j);
                x1 = loc2(k,1)-met(j);
                x2 = loc2(k,1)+met(j);
                sec2 = v2(y1:y2,x1:x2);
                score = abs(sec1-sec2);
                score = sum(sum(score));
                scrs(k,j) = score/((met(j)*2+1)^2);
                if (score<ms)
                    ms = score;
                    mind = k;
                    ttemp = score/((met(j)*2+1)^2);
                end
            end
            cor_loc(j,:) = [loc2(mind,:)-[met_m met_m]];
        else
            loc1(j,:) = [];
            cor_loc(j,:) = [];
            rem1 = rem1 +1;
        end
    end
    scrs = sort(scrs);
%     
%     figure
%     imshow(v1);hold on
%     plot(loc1(:,1),loc1(:,2),'r*','MarkerSize', 12); hold on
% %     plot(1500,481,'b*','MarkerSize', 12);
%     figure
%     imshow(vv2);hold on
%     plot(cor_loc(:,1),cor_loc(:,2),'r*','MarkerSize', 12); hold on
    newnp = np-rem-rem1;
    hmg_1 = ones(newnp,1);
    
    %%Homogenising and tranforming using hartleys normalization
    %%%% Frame 1
    l_p = [loc1 hmg_1];
    mu_x1 = sum(loc1(:,1))/newnp;
    mu_y1 = sum(loc1(:,2))/newnp;
    loc1(:,1) = loc1(:,1) - mu_x1;
    loc1(:,2) = loc1(:,2) - mu_y1;
    dummy1 = loc1.*loc1;
    dummy2 = sqrt(sum(dummy1,2));
    scale1 = sum(dummy2)/(newnp*sqrt(2));
    T_1 = [1/scale1 0 0;0 1/scale1 0; -mu_x1/scale1 -mu_y1/scale1 1];
%     T_1 = [1 0 0;0 1 0; -mu_x1 -mu_y1 1];
    loc1 = loc1/scale1;
%     loc1 = [loc1 hmg_1];
%     loc1 = loc1*T_1;
    

    %%%% Frame 2
    ch_p = [cor_loc hmg_1];
    mu_x2 = sum(cor_loc(:,1))/newnp;
    mu_y2 = sum(cor_loc(:,2))/newnp;
    cor_loc(:,1) = cor_loc(:,1) - mu_x2;
    cor_loc(:,2) = cor_loc(:,2) - mu_y2;
    dummy1 = cor_loc.*cor_loc;
    dummy2 = sqrt(sum(dummy1,2));
    scale2 = sum(dummy2)/(newnp*sqrt(2));
    T_2 = [1/scale2 0 0;0 1/scale2 0; -mu_x2/scale2 -mu_y2/scale2 1];
%     T_2 = [1 0 0;0 1 0; -mu_x2 -mu_y2 1];
    cor_loc = cor_loc/scale2;
%     ch_p = [cor_loc hmg_1]*T_2;
%     cor_loc = cor_loc*T_2;
    
    [fRANSAC, inl] = estimateFundamentalMatrix(loc1,cor_loc,'Method','RANSAC','NumTrials',2000,'DistanceThreshold',1e-4);    
    Ft(:,:,i-1) = fRANSAC;
    
    %%%Tranforming back the F matrix to scale up to remove normalization
    Ft(:,:,i-1) = T_1*Ft(:,:,i-1)*T_2';
    E = K'*Ft(:,:,i-1)*K;   % Converting to essential matrix
    
    [U,S,V] = svd(E);
    S(1,1) = 1;
    S(2,2) = 1;
    S(3,3) = 0;
    E = U*S*V';
    
    %% Getting the rotation and translation vectors from essential matrix
    W = [0 -1 0;1 0 0;0 0 1];
    if det(V) == -1
        V(:,3) = -V(:,3);
    end
    if det(U) == -1
        U(:,3) = -U(:,3);
    end
    t1 = U(:,3);
    t2 = -U(:,3);
    R1 = U*W*V';
    R2 = U*W'*V';
    P1 = [R1 t1];
    P2 = [R2 t1];
    P3 = [R1 t2];
    P4 = [R2 t2];
    sc = zeros(4,cfactor);
    
    %% validating the obtained projection matrices(needs to be corrected)
    for i = 1:cfactor%newnp
        t_p = ch_p(i,:)';
        t_p1 = l_p(i,:)';
        temp1 = pinv(K*P1)*t_p;
        temp0 = pinv(K*[eye(3) [0;0;0]])*t_p1;
        A = [temp0(1:3) -temp1(1:3)];
        B = P1(1:3,4);
        ot1 = pinv(A)*B;
        o1a = ot1(1)*temp0(1:4);
        o1b = K*P1*o1a;
        sc1 = o1a(3)*o1b(3)*min(o1a(3),o1b(3));
    %     sc(1,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
        temp2 = pinv(K*P2)*t_p;
        A = [temp0(1:3) -temp2(1:3)];
        B = P2(1:3,4);
        ot2 = pinv(A)*B;
        o2a = ot2(1)*temp0(1:4);
        o2b = K*P2*o2a;
        sc2 = o2a(3)*o2b(3)*min(o2a(3),o2b(3));
    %     sc(2,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
        temp3 = pinv(K*P3)*t_p;
        A = [temp0(1:3) -temp3(1:3)];
        B = P3(1:3,4);
        ot3 = pinv(A)*B;
        o3a = ot3(1)*temp0(1:4);
        o3b = K*P3*o3a;
        sc3 = o3a(3)*o3b(3)*min(o3a(3),o3b(3));
    %     sc(3,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
        temp4 = pinv(K*P4)*t_p;
        A = [temp0(1:3) -temp4(1:3)];
        B = P4(1:3,4);
        ot4 = pinv(A)*B;
        o4a = ot4(1)*temp0(1:4);
        o4b = K*P4*o4a;
        sc4 = o4a(3)*o4b(3)*min(o4a(3),o4b(3));
    %     sc(4,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
        tt = [sc1 sc2 sc3 sc4];
        [val,ind] = min(tt);
        sc(ind,i) = 1;
    end
    sc = sum(sc,2);
    
end
toc