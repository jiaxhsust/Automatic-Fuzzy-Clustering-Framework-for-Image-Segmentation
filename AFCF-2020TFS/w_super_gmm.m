function [Label,centers, U, obj_fcn,iter_n,Nan,covars,pi_k] =w_super_gmm(Super_L,Data, Number,cluster_n,options)
%% check  input variants
% if nargin ~= 2
%     error('Too many or too few input arguments!');
% end
if min(size(Data))<1
    error('Imput data of image, not an image ');
end
[data_n,dim]=size(Data);
GMM_WIDTH=1; 
Number=Number';
M_num=ones(cluster_n,1)*Number;  
default_options = [2;	% exponent for the partition matrix U
    100;	% max. number of iteration
    1e-5;	% min. amount of improvement
    0];	    % info display during iteration
if nargin ==4
    options = default_options;
else
    % If "options" is not fully specified, pad it with default values.
    if length(options) < 4
        tmp = default_options;
        tmp(1:length(options)) = options;
        options = tmp;
    end
    % If some entries of "options" are nan's, replace them with defaults.
    nan_index = find(isnan(options)==1);
    options(nan_index) = default_options(nan_index);
    if options(1) <= 1
        error('The exponent should be greater than 1!');
    end
end
expo = options(1);          % Exponent for U
max_iter = options(2);      % Max. iteration
min_impro = options(3);     % Min. improvement
display = options(4);       % Display info or not
obj_fcn = zeros(max_iter, 1);   % Array for objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nan=1;
Int=1;
while Nan
    [centers,r_ik] = FastFCM(Data,cluster_n,Number);  
    mf=M_num.*r_ik;
    for i = 1:cluster_n
        diffs = Data-(ones(data_n, 1) * centers(i,:));
        diffs = diffs.*(sqrt(mf(i,:))'*ones(1, dim));
        covars(:,:,i)=(diffs'*diffs)/sum(mf(i,:));              
        try
            if rank(covars(:,:,i)) < dim
                covars(:,:,i) = covars(:,:,i) + GMM_WIDTH.*eye(dim);
            end
        catch
            covars(:,:,i) = covars(:,:,i) + GMM_WIDTH.*eye(dim);
        end
    end
    
    A=sum(covars,3);
    Nan=sum(isnan(A(:)));
    Int=1+Int;
    if Nan>=1
        disp('There is NaN in covariance')
    end
    if Int==5;break;end   
end
if Nan~=0
    centers=[];
    U=[];
    obj_fcn=[];
    iter_n=[];
else
    pi_k=sum(mf,2)./sum(Number);
    MIN_COVAR = 1e-10;    % Minimum singular value of covariance matrix
    init_covars =covars;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for  w= 1:max_iter
        for i = 1:cluster_n
            diffs = Data-(ones(data_n, 1) * centers(i,:));
            mahalanobis2 = sum((diffs*inv(covars(:,:,i))) .* diffs,2);  
            dist(i,:) =pi_k(i).* (1/((2*pi)^(dim/2)*(det(covars(:,:,i)))^(1/2))).*exp(-mahalanobis2./2);  %
        end
        r_ik = dist./(ones(cluster_n,1)*sum(dist,1));  
        mf=M_num.*r_ik;
        centers = mf*Data./(sum(mf,2)*ones(1,size(Data,2))); %new center
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:cluster_n
            diffs = Data-(ones(data_n, 1) * centers(i,:));
            diffs = diffs.*(sqrt(mf(i,:))'*ones(1, dim));
            covars(:,:,i)=(diffs'*diffs)/sum(mf(i,:));
        end
        
        pi_k=sum(mf,2)./sum(Number);
        % Ensure that no covariance is too small
        A=sum(covars,3);
        Nan=sum(isnan(A(:)));
        if isreal(covars)==0 ||Nan,break; end
        
        for j = 1:cluster_n
            if min(svd(covars(:,:,j))) < MIN_COVAR
                covars(:,:,j) = init_covars(:,:,j);
            end
        end
        obj_fcn(w)= sum(log(sum(dist)));
        Uc{w}=r_ik;
        if w> 1
            if abs(max(max(Uc{w} - Uc{w-1}))) < min_impro, break; end
        end
    end
    iter_n=w;
    U=r_ik;
    obj_fcn=obj_fcn(1:w);
end


[~, GMM_L] = max(U, [], 1);
[rows,cols]=size(Super_L);
Label=zeros(rows,cols); % the output label
for i=1:max(Super_L(:))
    Label=Label+(Super_L==i)*GMM_L(i);
end
