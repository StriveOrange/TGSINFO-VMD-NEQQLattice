clc;
clear;
tent=0.5;
dim=1;
N=100;
TentX=rand(N,1);
for i=1:N
    for j=2:dim
        if TentX(i,j-1)<tent
            TentX(i,j)=TentX(i,j-1)/tent;
        elseif TentX(i,j-1)>=tent
            TentX(i,j)=(1-TentX(i,j-1))/(1-tent);
        end
    end
end
plot(TentX,'.')
xlabel('维度')
ylabel('混沌值')
figure
hist(TentX)
xlabel('混沌值')
ylabel('频数')