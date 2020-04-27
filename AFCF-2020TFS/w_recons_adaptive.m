%% MMGR-WT achieves superpixel segmentation using multiscale morphological gradient reconstruction
% L_seg is lable image with line
% i is the maximal iterations
% diff is the difference between the previous result and current gradient
function [L_seg,f_g2,diff,i]=w_recons_adaptive(f,se_start, options) 
% sto_con denotes end condition
%       OPTIONS(1): maximum number of iterations          (default: 15)
%       OPTIONS(3): minimum amount of improvement         (default: 1e-4)
%sto_con=[15,0.0001];
if nargin ~= 2 & nargin ~= 3
	error('Too many or too few input arguments!');
end
% Change the following to set default options
default_options = [8;	% max. number of iteration
		1e-1];	% min. amount of improvement
if nargin == 2
	options = default_options;
    max_itr= default_options(1);
	min_impro= default_options(2);
else
	% If "options" is not fully specified, pad it with default values.
		max_itr= options(1);
		min_impro=options(2);
end
% max_itr=sto_con(1); 
% min_impro=sto_con(2);
%% multi-scale morphological gradient reconstruction
f_g=zeros(size(f,1),size(f,2));diff=zeros(max_itr,1);
for i=se_start:max_itr
    gx=w_recons_CO(f,strel('disk',i+se_start-1)); %gradient reconstruction with different SE
    f_g2=max(f_g,double(gx));%max corresponding to zeros in f_g; min corresponding to ones in f_g
    f_g1=f_g;f_g=f_g2;
% diff(i)=abs(double(max(max(watershed(f_g1))))-double(max(max(watershed(f_g2)))));
% check termination condition,2;
    diff(i)=max(max((abs(f_g1 - f_g2))));
	if i > 1
		if diff(i) < min_impro, break; end
    end  
    %diff(i)
end
%% step4  watershed
L_seg=watershed(f_g2);
%MMR_seg=Label_image_WT(f,L_seg);