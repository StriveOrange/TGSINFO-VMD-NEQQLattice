%% Population initialization process improved based on Tent chaos mapping.
function X=initialization(nP,dim,ub,lb)
TentX = rand(nP,dim);
tent=1-eps;
for i=1:nP
    for j=2:dim
        if TentX(i,j-1)<tent
            TentX(i,j)=TentX(i,j-1)/tent;
        elseif TentX(i,j-1)>=tent
            TentX(i,j)=(1-TentX(i,j-1))/(1-tent);
        end
    end
end
X = lb + TentX.*(ub - lb);
end