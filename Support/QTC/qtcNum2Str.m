function qtcStr=qtcNum2Str(qtcNum)
qtcStr=cell(size(qtcNum,1),1);
for (i=1:size(qtcNum,1))
    for (j=1:size(qtcNum,2))
        s(j)='0';
        switch qtcNum(i,j)
            case 0
                s(j)='0';
            case 1
                s(j)='+';
            case -1
                s(j)='-';
        end;
    end;
    qtcStr{i}=s;
end;