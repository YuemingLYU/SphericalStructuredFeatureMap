function [ Y ] = SSF( X,n,base )
[d,N]= size(X);
X= X(1:d/2,:) +1i*X(d/2+1:end,:);
S= (zeros(n/2,N));
S(1+base,:)=X;
clear X;
S=fft(S);
Y =  1/sqrt(d/2)*S;
clear S
Y= [real(Y);imag(Y)];

end

