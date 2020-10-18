function [base,bestD] = LogEnergyOP_Demo(s, n,Tn)
str= strcat(['./results/LogEnergyOP_Demo_',num2str(s),'_n_', num2str(n),'_T_', num2str(Tn)]);
try
     load(str);
catch

tic;

         hist=zeros(1,Tn);

    base = randperm(n-1,s);
        
     
      T=mod((1:n-1)'*base/n,1);
     Te=cos(2*pi*T)+ 1i*sin(2*pi*T);
     z= sum(Te,2);  
     
 bestD=-(sum(log(((1-1/(s^2)*imag(z).^2)))) + sum(log(((1-1/(s^2)*real(z).^2)))));
  
    B= round(5*10^7/n);
    TB = floor(n/s);
    TN = 5;
  
    clear T Te;
    
   
    
    for k=1:Tn
       % tic
        if mod(k,1)==0
           disp(strcat([num2str(k),' bestD: ',num2str(bestD,10)]));
            hist(k)= bestD;
        end

        
        for ii=1:s
   
             bestD=-(sum(log(((1-1/(s^2)*imag(z).^2)))) + sum(log(((1-1/(s^2)*real(z).^2))))) ;

             Ti=mod((1:n-1)'*base(ii)/n,1);
             Ti=cos(2*pi*Ti) + 1i*sin(2*pi*Ti);
            t= z - Ti;


            
            
              Block = randperm((n-1),min(B,(n-1)));

            Block = setdiff(Block,base);

             tmp = mod((1:n-1)'*Block/n,1);
                       

             Ti= cos(2*pi*tmp) + 1i*sin(2*pi*tmp);
             Dis= bsxfun(@plus,Ti,t);
                 
             Mdis = -  ( sum(log((1-1/(s^2)*real(Dis).^2)),1) + sum(log((1-1/(s^2)*imag(Dis).^2)),1) ) ;
             [Tva,Tid]=min(Mdis);


                if Tva<bestD
                   bestD=Tva;
                   base(ii)=Block(Tid);
                  
              
                end    


            Ti=mod((1:n-1)'*base(ii)/n,1);
            Ti= cos(2*pi*Ti) + 1i*sin(2*pi*Ti);

            z= t + Ti;

        end
  
    end

   disp(strcat(['bestDistance: ',num2str(bestD,10)]));


clear Ti  z t tmp Dis



  save(str,'hist','bestD','base');
          
          
    
          
toc;
end

load(str);



end



