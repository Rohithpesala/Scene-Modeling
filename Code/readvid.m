v = VideoReader('VID_20151012_141216.mp4');
nf = v.NumberOfFrames;
h = v.Height;
w = v.Width;
np = 25;            %%% number of points to use for correspondences
Ft = zeros(3,3,nf-1);
tfactor = 2;      %%% tolerance factor

%% Calibrated intrinsic params
K = [2088 0 1269; 0 2066 969;0 0 1];


%%%% Reading initial frame %%%%
% v1 = rgb2gray(read(v,1));
% points = detectSURFFeatures(v1);
% p = points.selectStrongest(1);
% p(1).Location
% imshow(v1); hold on;
% plot(p);
tic
%%% looping over each frame
for i = 32:32
    i
    %%% two frames from video
    v1 = im2double(rgb2gray(read(v,i-1)));
    v2 = im2double(rgb2gray(read(v,i+5)));
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
    p1 = points1.selectStrongest(np);
    loc1 = zeros(np,2);
    met = zeros(np,1);
    scrs = zeros(np*3,np);
    cor_loc = zeros(np,2);
    for j = 1:np
        loc1(j,:) = round(p1(j).Location);
        met(j) = round(p1(j).Scale);
    end
    met = met*6;
    
    %% adding padding for second frame with max radius to avoid conflicts. Please neglect this padding
    met_m = ceil(max(met));
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
    for j = 1:np
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
    end
    scrs = sort(scrs);
    rem = 0;
    for j = 1:np
        if scrs(2,j) <= scrs(1,j)*tfactor
            loc1(j-rem,:) = [];
            cor_loc(j-rem,:) = [];
            rem = rem+1;
        end
    end
    figure
    imshow(v1);hold on
    plot(loc1(:,1),loc1(:,2),'r*','MarkerSize', 12);
    figure
    imshow(vv2);hold on
    plot(cor_loc(:,1),cor_loc(:,2),'r*','MarkerSize', 12);
%     cor_loc = loc1;  %%% remove after experimentation
%     cor_loc(:,2) = cor_loc(:,2)+30;
    %%%% Transforming FPs %%%%
%     loc1 = po1';
%     cor_loc = po2';
    newnp = np-rem;
    hmg_1 = ones(newnp,1);
    
    %%Homogenising and tranforming using hartleys normalization
%     tloc1 = [loc1 hmg_1];
    %%%% Frame 1
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
    loc1 = [loc1 hmg_1];
%     loc1 = loc1*T_1;
    

    tloc1 = [cor_loc hmg_1];
    %%%% Frame 2
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
    cor_loc = [cor_loc hmg_1];
%     cor_loc = cor_loc*T_2;
    
    %%%% Method 2: searching full space for correspondences %%%%
%     for j = 1:np
%         sec1 = v1(loc(j,2)-met(j) : loc(j,2)-met(j), loc(j,1)-met(j) : loc(j,1)-met(j));
%         %sm = zeros(h,w);
%         min_sc = inf;
%         min_x = 1;
%         min_y = 1;
%         for k = met(j)+1:h-met(j)-1
%             k
%             for l = met(j)+1:w-met(j)-1
%                 sec2 = v2(k-met(j):k+met(j), l-met(j):l+met(j));
%                 score = abs(sec1-sec2);
%                 score = sum(sum(score));
%                 if (score<min_sc)
%                     min_sc = score;
%                     min_x = l;
%                     min_y = k;
%                 end
%             end                       
%         end
%         cor_loc(i,:) = [min_x min_y];
%     end

    %%% Forming the U matrix
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
    [U,S,V] = svd(F);
    S(3,3) = 0;         % equating the last singular value to 0 as F is rank 2
    Ft(:,:,i-1) = U*S*V';
    
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
    t1 = U(:,3);
    t2 = -U(:,3);
    R1 = U*W*V';
    R2 = U*W'*V';
    P1 = [R1 t1];
    P2 = [R2 t1];
    P3 = [R1 t2];
    P4 = [R2 t2];
    sc = zeros(4,newnp);
    
    %% validating the obtained projection matrices
    for i = 1:newnp
        t_p = cor_loc(i,:)';
        temp = pinv(P1)*t_p;
        sc(1,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
        temp = pinv(P2)*t_p;
        sc(2,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
        temp = pinv(P3)*t_p;
        sc(3,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
        temp = pinv(P4)*t_p;
        sc(4,i) = temp(4)*temp(3)/abs(temp(4)*temp(3));
    end
end
toc
