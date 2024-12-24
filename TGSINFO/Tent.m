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
xlabel('Dimensions')
ylabel('Chaotic values')
figure
hist(TentX)
xlabel('Chaotic values')
ylabel('Frequency')
