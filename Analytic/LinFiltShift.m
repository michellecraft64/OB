function [kFilt,Shft,Coeff,lsApprx]=LinFiltShift(n,st1,y,trnc)
% fit using least-squares; linear filter going back 'n' steps
% !! ASSUMING n is an int, st1 is a COLUMN vector!!

Lt=length(st1); %assuming all same size
idSpnEnd=floor(1/3*(Lt-1))+1;  % approx when spont ends, to fill-in 1:n-1

%-- OUTPUTS-- 
    Coeff=zeros(2,1);
    kFilt=zeros(n,1);
    Shft=0;
    lsApprx=zeros(Lt,1);
    
    if(n>Lt)
        disp('Need n to be less than length of vecs!');
        return;
    end
    sMatr=toeplitz(st1(n:Lt),st1(n:-1:1)');
    
    [Qm,Rm]=qr([sMatr ones(Lt-n+1,1)],0);
    lsc=inv(Rm)*Qm'*y(n:Lt);
    kFiltTrunc=lsc(1:round(n*trnc));
    kFiltAdd=zeros(n-length(kFiltTrunc),1);
    kFilt=[kFiltTrunc;kFiltAdd]; Shft=lsc(end);
    
    kMat=sparse(convmtx(kFilt(end:-1:1)',Lt-n+1));
    lsApprx(n:Lt)=kMat*st1+Shft; %perform convolution and shift
    if(n<idSpnEnd)
        lsApprx(1:n-1)=mean(lsApprx(n:idSpnEnd)); %take avg up to spont 
    else
        lsApprx(1:n)=lsApprx(n+1); 
    end
end