function [c,n]=qtc2case(qtc);
dimension=size(qtc,2);
mult=3.^((dimension-1):-1:0);
coded=(qtc+1).*repmat(mult,size(qtc,1),1);
n=(sum(coded,2))'+1;
c=char(n+32);