function [ SSFmaxError,SSFmeanError ,RFFmaxError,RFFmeanError] = Test(n,base,dim ,M, name, sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


 N=2000;
 
load(name,'Xdata');
[NN,d] = size((Xdata));
 Xdata=Xdata/sigma;


Norm1 = sqrt(icdf('chi2',1/2,d) );  % For Gaussian kernel 
Norm2= sqrt(d);                     % For cosKerel and angelKernel from the closed form

    
    
for ii=1:M
  
%  
ii

id = randperm(NN,N);
X=full(Xdata(id,:))';

  d= size(X,1);
 D0 = sign(randn(d,1));

 X =bsxfun(@times,X,D0);

Nn=n;

% Construct SSF
 [ZXn] = SSF( X,n,base );
 

 
 %% Gaussian Kernel
 
 % SSF
 ZX = Norm1 *ZXn;
 Y=[cos(ZX);sin(ZX)] * sqrt(1/Nn);
 T=Y'*Y;
 
 
 % RFF
  B= randn(Nn,dim) ;
  BX=B*X;
 Y2=[cos(BX);sin(BX)] * sqrt(1/Nn);
 T2=Y2'*Y2;
 
 
 K=GaussianKernl(X);
 
 
 TK1= abs(T-K);
 TK2=abs(T2-K);
 
 a=max(TK1(:));
 b=max(TK2(:));

t= max(K(:));
SSFmaxError(ii,1)=a/t;
 RFFmaxError(ii,1)=b/t;
 
disp(strcat(['GaussianKernelMaxErrorSSF: ',num2str(SSFmaxError(ii,1),10),' RFF: ',num2str(RFFmaxError(ii,1),10)],' RDiff: ', num2str((a-b)/b,10)));


t= sqrt(sum(K(:).^2));
  a= sqrt(sum(TK1(:).^2));
  b= sqrt(sum(TK2(:).^2));

 SSFmeanError(ii,1)=a/t;
 RFFmeanError(ii,1)=b/t;
 

 disp(strcat(['GaussianKernelMeanErrorSSF: ',num2str(SSFmeanError(ii,1),10),' RFF: ',num2str(RFFmeanError(ii,1),10)],' RDiff: ', num2str((a-b)/b,10)));
 
 
 
 %% cosKerel
 
 
% SSF
ZX = Norm2 * ZXn;
 Y= [(ZX>0);(ZX<=0)] * sqrt(1/Nn);
 T=Y'*Y;
 
 %RFF
 Y2=[(BX>0);(BX<=0)] * sqrt(1/Nn);
 T2=Y2'*Y2;
 
 
 
 K=cosKerenl(X);
 
 
 TK1= abs(T-K);
 TK2=abs(T2-K);
 
 a=max(TK1(:));
 b=max(TK2(:));
 
 
t= max(K(:));
SSFmaxError(ii,2)=a/t;
RFFmaxError(ii,2)=b/t;

disp(strcat(['CosKernelMaxErrorSSF: ',num2str(SSFmaxError(ii,2),10),' RFF: ',num2str(RFFmaxError(ii,2),10)],' RDiff: ', num2str((a-b)/b,10)));

 t= sqrt(sum(K(:).^2));
  a= sqrt(sum(TK1(:).^2));
  b= sqrt(sum(TK2(:).^2));

 SSFmeanError(ii,2)=a/t;
 RFFmeanError(ii,2)=b/t;
 
 disp(strcat(['CosKernelMeanErrorSSF: ',num2str(SSFmeanError(ii,2),10),' RFF: ',num2str(RFFmeanError(ii,2),10)],' RDiff: ', num2str((a-b)/b,10)));
 
 
 
 
 %% angelKernel
 
% SSF
Y= [max(ZX,0);max(-ZX,0)] * sqrt(1/Nn);
 T=Y'*Y;
 
 % RFF
 Y2=[max(BX,0);max(-BX,0)] * sqrt(1/Nn);
 T2=Y2'*Y2;
 
 

 
 K=angleKernel(X);
 
 
 TK1= abs(T-K);
 TK2=abs(T2-K);
 
 a=max(TK1(:));
 b=max(TK2(:));
 
t= max(K(:));
SSFmaxError(ii,3)=a/t;
 RFFmaxError(ii,3)=b/t;
 
disp(strcat(['AngleKernelMaxErrorSSF: ',num2str(SSFmaxError(ii,3),10),' RFF: ',num2str(RFFmaxError(ii,3),10)],' RDiff: ', num2str((a-b)/b,10)));

t= sqrt(sum(K(:).^2));
  a= sqrt(sum(TK1(:).^2));
  b= sqrt(sum(TK2(:).^2));
  

 SSFmeanError(ii,3)=a/t;
 RFFmeanError(ii,3)=b/t;

disp(strcat(['AngleKernelMeanErrorSSF: ',num2str( SSFmeanError(ii,3),10),' RFF: ',num2str(RFFmeanError(ii,3),10)],' RDiff: ', num2str((a-b)/b,10)));
 
 
 

 end




end

function [K]= GaussianKernl(X)

 D=pdist2(X',X','euclidean');
 
 K = exp(-D.^2/2);

end


function [K]=cosKerenl(X)

  %  A=1-pdist2(X',X','cosine');
  normX = sqrt(sum(X.^2,1));
  X= bsxfun(@rdivide,X,normX);
  A= X'*X;
  A= min(A,1);
  A=max(A,-1);
    K= 1- acos(A)/pi;

end


function [K] = angleKernel(X)

normX = sqrt(sum(X.^2,1));
X= bsxfun(@rdivide,X,normX);
 A= X'*X;
 A= min(A,1);
A=max(A,-1);
theta = acos( A);

T= sin(theta) + (pi-theta).*cos(theta);
K = bsxfun(@times, bsxfun(@times,T,normX), normX')/pi;

end



