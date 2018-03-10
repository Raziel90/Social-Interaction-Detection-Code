function hmm = qtcCND;

qtc=[];
for i1=1:3
	for i2=1:3
		for i3=1:3
			for i4=1:3
				qtc(end+1,:)=[i1-2 i2-2 i3-2 i4-2];
			end;
		end;
	end;
end;

% Find valid transitions
% Valid transitions are represented by a 1 in d; 0 or 2 are invalid
validqtc=[];
d=zeros(83,83);
for i1=1:size(qtc,1)
    for i2=1:size(qtc,1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % self transition = 0, transition form - to + and vice versa = 2
        % transitions from - to 0 or + to 0 and vice versa = 1
        d(i1+1,i2+1)=max(abs(qtc(i1,:)-qtc(i2,:)));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find invalid transitions according to CND:
        % 1,2: -000 <> 0-00 | +000 <> 0+00 | -000 <> 0+00 | +000 <> 0-00
        % 1,3: -000 <> 00-0 | +000 <> 00+0 | -000 <> 00+0 | +000 <> 00-0
        % 1,4: -000 <> 000- | +000 <> 000+ | -000 <> 000+ | +000 <> 000-
        % 2,3: 0-00 <> 00-0 | 0+00 <> 00+0 | 0-00 <> 00+0 | 0+00 <> 00-0
        % 2,4: 0-00 <> 000- | 0+00 <> 000+ | 0-00 <> 000+ | 0+00 <> 000-
        % 3,4: 00-0 <> 000- | 00+0 <> 000+ | 00-0 <> 000+ | 00+0 <> 000-
        if d(i1+1,i2+1)==1
            for j1=1:size(qtc(i1,:),2)-1
                for j2=j1+1:size(qtc(i2,:),2)
                    if sum(abs(qtc(i1,[j1,j2])))==1 ...
                            & sum(abs(qtc(i2,[j1,j2])))==1
                        if max(abs(qtc(i1,[j1,j2])-qtc(i2,[j1,j2]))) > 0 ...
                                & sum(qtc(i1,[j1,j2])-qtc(i2,[j1,j2]))~=1
                            d(i1+1,i2+1)=2;
                            break;
                        end
                    end
                end
                if d(i1+1,i2+1) ~= 1
                    break;
                end
            end
        end
        % Create list of valid transitions
        if d(i1+1,i2+1)==1
            validqtc(end+1,:)=[qtc(i1,:), 5 ,qtc(i2,:)];
        end
    end
end

validtr=double(d==1);

validtr(1,:)=1;
validtr(:,1)=0;
validtr(:,end)=1;
validtr(1,end)=0;
validtr(end,:)=0;

validtr=validtr+eye(83,83)*0.00001;
validtr(1,1)=0;
t=validtr./repmat(sum(validtr,2),1,size(validtr,2));

%emmissions=ones(83,83)*(0.001/82)+eye(83,83)*0.999;
emmissions=eye(83,83);

% emmissions(1,2:end)=0;
%emmissions(end,1:end-1)=0.0001;
% emmissions=emmissions./repmat(sum(emmissions,2),1,size(emmissions,2));

hmm.t=t;
hmm.e=emmissions;
hmm.q = validqtc;