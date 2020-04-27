function  [center,U,obj_fcn]=FastFCM(data,cluster,Number)
if nargin ~= 3
    error('Too many or too few input arguments!');
end
% Change the following to set default options
options = [2;   % exponent for the partition matrix
    100;    % max. number of iteration
    1e-5;   % min. amount of improvement
    1]; % info display during iteration
expo = options(1);      % Exponent for U
max_iter = options(2);      % Max. iteration
min_impro = options(3);     % Min. improvement
display = options(4);       % Display info or not
obj_fcn = zeros(max_iter, 1);   % Array for objective function
M_num=ones(cluster,1)*Number;
U=initfcm(cluster,length(data));   % Initial fuzzy partition
for w= 1:max_iter
    mf=M_num.*(U.^expo);
    center = mf*data./((ones(size(data, 2), 1)*sum(mf'))');
    dist = distfcm(center, data);       % fill the distance matrix
    tmp = dist.^(-2/(expo-1));      % calculate new U, suppose expo != 1
    U= tmp./(ones(cluster, 1)*sum(tmp));
    obj_fcn(w)=sum(sum((dist.^2).*mf));
    %     if display,
    %          fprintf('Iteration count = %d,obj_fcn = %f\n',w,obj_fcn(w));
    %      end
    if w > 1
        if abs(obj_fcn(w)-obj_fcn(w-1)) < min_impro, break; end
    end
end



