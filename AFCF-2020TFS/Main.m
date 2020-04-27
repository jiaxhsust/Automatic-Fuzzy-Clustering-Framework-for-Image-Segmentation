% Please kindly cite the paper Tao Lei, Peng Liu, Xiaohong Jia, Xuande Zhang, Hongying Meng and Asoke K. Nandi, 
%"Automatic Fuzzy Clustering Framework for Image Segmentation," 
%IEEE Transactions on Fuzzy Systems,2019,Doi:10.1109/TFUZZ.2019.2930030

% The code was written by Tao Lei, Peng Liu and Xiaohong Jia in 2019.
%%%% Welcome to our Research Group website:https://aimv.sust.edu.cn/lwcg.htm
clc
clear all
close all
f_ori=imread('.\Images\100007.jpg'); 
figure,imshow(f_ori);
[rows,cols,dim]=size(f_ori);
%% parameters computation
f_size=min(size(f_ori,1),size(f_ori,2));
gauF_w=2*ceil(f_size/500)+1;rs=ceil(f_size/500)+1;rm=ceil(f_size/100)+5;
%% gaussian filtering
sigma=1.0;gausFilter=fspecial('gaussian',[gauF_w gauF_w],sigma);g=imfilter(f_ori,gausFilter,'replicate');
%% Actually, you can use SE (2015-PAMI) to obtain better gradient images
F_ori=rgb2lab(f_ori); 
a1=sgrad_edge(normalized(F_ori(:,:,1))).^2;b1=sgrad_edge(abs(normalized(F_ori(:,:,2)))).^2;c1=sgrad_edge(normalized(F_ori(:,:,3))).^2;
Gradient=sqrt(a1+b1+c1); 
%% AMR-WT, you can find detailed technology and algorithm in TIP-2019
%Please cite the paper "Tao Lei, Xiaohong Jia,Tongliang Liu,Shigang Liu,Hongying Meng,and Asoke K. Nandi, 
%Adaptive Morphological Reconstruction for Seeded Image Segmentation,
%IEEE Transactions on Image Processing, vol.28, no.11, pp.5510-5523, Nov. 2019."


Super_L=w_recons_adaptive(Gradient,rs,[rm 0.0001]); 
L2=imdilate(Super_L,strel('square',2));
[Lseg,~,Num,centerLab]=Label_image(f_ori,L2);
BW = boundarymask(L2);
Super_line = imoverlay(Lseg,BW,'red');figure,imshow(Super_line);
%% Desnity peaks--2014Science
percent=2;
[Lab2,gamma,rho,delta,cluster_n,icl,Global_C_Lab,Cluster_Number]=w_DPRS_interval(centerLab,Num,percent);
%% Gaussian Mixed Model
[final_Label,center_Lab] =w_super_gmm(L2,centerLab,Num,cluster_n);
fs=Label_image(f_ori,final_Label);
BW = boundarymask(final_Label);
Lseg_line = imoverlay(fs,BW,'red');figure,imshow(Lseg_line);

%% drawing figures 
figure,plot(gamma,'s')
hold on
plot(gamma(1:cluster_n),'s','MarkerFaceColor','r','Linewidth',2.5,'MarkerEdgeColor','r')  % clustering centers are 
xlabel ('n');
ylabel ('\phi');
