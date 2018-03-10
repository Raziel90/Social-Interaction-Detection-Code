function label=qtcCase2Label(c)
for k=1:length(c)
    q=[27 9 3 1];

    rc=c(k)-1;

    f=floor(rc(1)/q(1));
    r=rem(rc(1),q(1));

    for i=2:length(q)
        rc(i)=rc(i-1)-f(i-1)*q(i-1);
        f(i)=floor(rc(i)/q(i));
        r(i)=rem(rc(i),q(i));
    end;
    qtc(k,:)=f-1;
   
end
label=qtcNum2Str(qtc);

