v = VideoReader('VID_20151012_141216.mp4');
nf = v.NumberOfFrames;
h = v.Height;
w = v.Width;
np = 10;            %%% number of points to use for correspondences
%%%% Reading initial frame %%%%
v1 = rgb2gray(read(v,1));
% p(1).Location
% imshow(v1); hold on;
% plot(p);
tic
for i = 2:2
    v2 = rgb2gray(read(v,i));
    points = detectSURFFeatures(v1);
    p = points.selectStrongest(np);
    loc = zeros(np,2);
    met = zeros(np,1);
    cor_loc = zeros(np,2);
    for j = 1:np
        loc(j,:) = round(p(j).Location);
        met(j) = round(p(j).Scale);
    end
    for j = 1:np
        sec1 = v1(loc(j,2)-met(j) : loc(j,2)-met(j), loc(j,1)-met(j) : loc(j,1)-met(j));
        %sm = zeros(h,w);
        min_sc = inf;
        min_x = 1;
        min_y = 1;
    end
    for k = met(j)+1:h-met(j)-1
        k
        for l = met(j)+1:w-met(j)-1
            for j = 1:np
                
            sec2 = v2(k-met(j):k+met(j), l-met(j):l+met(j));
            score = abs(sec1-sec2);
            score = sum(sum(score));
            if (score<min_sc)
                min_sc = score;
                min_x = l;
                min_y = k;
            end
            end
        end                       
    end
    cor_loc(i,:) = [min_x min_y];
    
    A = zeros(np,9);
    
end
toc

