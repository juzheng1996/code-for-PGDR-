clear
clc

%Data loading
n=52;
load('taoyuan2023.mat');
for i=1:52
    dataset(i,:)=sum(taoyuan2023((i-1)*7+1:i*7,:));
end

%Variable assignment
y1 = dataset(:,1);
y2 = dataset(:,2:3);
x = dataset(:,4:end);

%calculate
[directions, pgdr_variable] = pgdr(y1,y2,x);







