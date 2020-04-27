function [Lab2,gamma2,rho,delta,NCLUST,icl,Global_C,Cluster_Number]=w_DPRS_interval(f,Number,percent)
%% parameters
% f is an input matrix of size M*N
% the default value of percent is 2
% Lab2 is the classification Labels
% gamma is the decision graph
% rho is local density
% delta the distance of local density
% NCLUST is the number of clusters that we estimate
% icl is the clustering center pos
% addpath('.\Tools\')
ND=size(f,1); % the number of samples
Dist=pdist(double(f));   % the Euclidean distance between pairs of objects in M*N matrix
dist=double(squareform(Dist)); % together with the function of pdist to obtain a squareform distcance
%figure,imshow(mat2gray(dist)) % showing the matrix of dist
N=size(Dist,2);
if N<=50
    N=50;
end
position=round(N*percent/100);  % computing the threshold used in dc, dc is robust to samples
sda=sort(Dist(:)); % ordering for Dist in order to compute delta comveniently
dc=sda(position).^2; %dc is a cutoff distance
if dc==0
    xx=sda(sda>0);
    dc=xx(1).^2;
end
%% computing local density rho
rho=zeros(1,ND);delta=rho;
for i=1:ND
    rho(i)=sum(exp(-dist(i,:).^2./dc).*Number')-1;
end

maxd=max(max(dist));
%% computing delta
[~,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;
for i2=2:ND
    delta(ordrho(i2))=maxd;
    for j2=1:i2-1
        if(dist(ordrho(i2),ordrho(j2))<delta(ordrho(i2))) % searching the minimal distance using loop programe
            delta(ordrho(i2))=dist(ordrho(i2),ordrho(j2));
            nneigh(ordrho(i2))=ordrho(j2); % record the nearest data which has a higher density peaks
        end
    end
end
delta(ordrho(1))=max(delta(:));
gamma=(rho).*(delta);
%% automatically compute k, if we implement kmeans on gamma, we will divide gamma into two classifications.
%the one is peaks, the other is non-peaks. The number of peaks is k
[gamma2,idx2]=sort(gamma,'descend');
Gamma=normalized(gamma2);
figure,plot(Gamma,'s')
xlabel ('n')
ylabel ('\gamma')

Data=Gamma';
psi=1000;   
eta=0.1;    
Ndata = Rescale(psi,eta,Data);  
num=size(Ndata,1);
if num>2
    zr = 10e-11;    
    [gamma2,idx2]=sort(Ndata','descend');
    ad=abs(diff(gamma2));  
    gamma2(gamma2<zr)=eps;
    ad(1)=0;  
    if num>20
        ad = ad(1:floor(0.9*end));
    end 
    [~, cs]=sort(ad,'descend'); %
    %     sprintf('Suggested cluster number is: %d, %d, %d, %d, %d', cs(1))
    NCLUST=cs(1);
else
    NCLUST=2;
end

if NCLUST==1
    NCLUST=2;
end
Cluster_Number=[NCLUST];
%% automatically choose clustering center
Lab2=-1*ones(1,ND); 
for i3=1:NCLUST
    Lab2(idx2(i3))=i3;
end
icl=idx2(1:NCLUST);
Global_C=f(icl,:);
%% classification, assignation
for i4=1:ND
    if (Lab2(ordrho(i4))==-1)
        if nneigh(ordrho(i4))==0
            Lab2(ordrho(i4))=1;
        else
            Lab2(ordrho(i4))=Lab2(nneigh(ordrho(i4)));
        end
    end
end
