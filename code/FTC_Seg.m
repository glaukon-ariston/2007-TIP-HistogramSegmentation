%%%%%%%%%%%%%%%%%
% The following functions implement the Fine to Coarse Histogram Segmentation described in 
% // J. Delon, A. Desolneux, J-L. Lisani and A-B. Petro, A non
% parametric approach for histogram segmentation, IEEE Transactions
% on Image Processing, vol.16, no 1, pp.253-261, Jan. 2007. //
%
% Usage :
% u = double(imread('../images/lena.png'));
% H = hist(u(:),0:255);
% idx=FTC_Seg(H,0);
% idx should contain the list of all minima separating the modes of
% H 

%   Copyright (c) 2016 Julie Delon
%%%%%%%%%%%%%%%%%

function idx = FTC_Seg(H,e)

% FTC_seg

% H = the integer-valued histogram to be segmented. 

% e = parameter of the segmentation 
% (corresponds to e = -log10(epsilon) in the paper)
% large e => coarse segmentation
% small e => fine segmentation


logId = fopen('val_matlab.m', 'w');
lH = length(H);

%% find the list of local minima and maxima of H
[p,idx_max] = findpeaks(H);
% https://www.mathworks.com/matlabcentral/answers/44227-finding-local-minimums-maximums-for-a-set-of-data
Hinv = 1.01*max(H) - H;
[p,idx_min] = findpeaks(Hinv);
idx = sort([idx_min,idx_max]);
if (idx(1)~=1) idx = [1,idx]; end
if (idx(end)~=lH) idx = [idx,lH]; end

%% from python/histogramSegmentation.py (findpeaks in Python get different results than in MatLab)
% idx = [1 48 49 50 67 68 69 71 72 73 74 76 77 79 81 82 83 84 85 92 93 100 102 103 110 112 113 114 116 119 120 127 128 130 136 137 138 145 147 148 150 155 169 170 171 173 176 178 186 188 189 196 197 199 200 202 203 206 207 209 211 212 221 222 234 235 237 239 241 243 244 245 256];

% find if idx starts with a minimum or a maximum
 if H(idx(1)) < H(idx(2))
    begins_with_min =  1;
 else
    begins_with_min =  0;
 end  

%% FILL THE LIST OF ENTROPIES FOR ALL MODES 
% The merging of two contiguous modes [a,b] and [b,c] can be done in two ways, 
% either by using the maximum M1 on [a,b] and by testing the decreasing hypothesis on [M1,c], 
% or by using the maximum M2 on [b,c] and by testing the increasing hypothesis on [a,M2]. 
% For each configuration, we compute the entropy of the worst interval against the considered hypothesis.

K = length(idx);
val=zeros(1,K-3);

% Loop on all optimas
for k =1:K-3
    
   % decide if we want to test the increasing or decreasing hypothesis on
   % [idx(k),idx(k+3)]
   
   if ((begins_with_min && mod(k,2)==1) || (~begins_with_min && mod(k,2)==0))
       inc = 1; 
   else inc = 0; 
   end
   
   % compute the max entropy on the interval [k,k+3]
  val(k) = max_entropy(H,idx(k),idx(k+3),e,inc);
      
end

%%%%%%%% MERGING of MODES

[valmin, kmin] = min(val); %[idx(kmin), idx(kmin+3)] is the first interval to merge
fprintf(logId, "val[%d] = [%s]\n", length(val), sprintf("%f ", val))

while(~isempty(val) && valmin<0)
    % update the list of min, max
    idx = [idx(1:kmin),idx(kmin+3:end)];       
    val = [val(1:min(kmin,end)),val(kmin+3:end)];
    val = val(1:length(idx)-3);
    fprintf(logId, "val[%d] = [%s]\n", length(val), sprintf("%f ", val))
    % update max_entropy around the removed optima 
    for j=max(kmin-2,1):min(kmin,length(val))
        % decide if increasing or decreasing
        if (begins_with_min && mod(j,2)==1) || (~begins_with_min && mod(j,2)==0)
            inc = 1; 
        else inc = 0; 
        end
        % update the max entropy on the interval [k,k+3]
        val(j) = max_entropy(H,idx(j),idx(j+3),e,inc);
        fprintf(logId, 'j %d inc %d val(j) %f\n', j, inc, val(j))
    end
   [valmin, kmin] = min(val);
end

if (begins_with_min) 
    idx = idx(1:2:end); 
else
    idx = idx(2:2:end);     
end
fprintf(logId, "idx[%d] = [%s]\n", length(idx), sprintf("%d ", idx))
fclose(logId);

%% Display the segmentation
bar(H,':r');
hold on;
for k = 1:length(idx)
    line([idx(k) idx(k)], [0 max(H(:))]);
end
hold off;


end





