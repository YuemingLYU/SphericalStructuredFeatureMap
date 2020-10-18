


 name='mnist';dim=784; sigma=15;
 n= 2*dim ; % number of feature maps. A pair [cos, sin] is a feature

 M=30   % number of independent runs
 Tn=5;
 s=(n-1)/2;

  
%%

% searching base vector for SSF
[base,bestD] = LogEnergyOP_Demo(dim/2, ceil(n/2),Tn);


%%
 
 
 [ SSFmaxError,SSFmeanError ,RFFmaxError,RFFmeanError] = Test( n,base,dim,M, name,sigma );


 
 %% approximation
 
  SSFmax = mean((SSFmaxError),1)
 
 SSFMean= mean((SSFmeanError),1)
 
 
 RFFmax = mean((RFFmaxError),1)
 
 RFFMean= mean((RFFmeanError),1)
 

 
 
 str = strcat(['./err_',name,'_dim_', num2str(dim), '_n_',num2str(n/dim),'_M_',num2str(M)]);
 
 save(str,'SSFmax','SSFMean','RFFmax','RFFMean');