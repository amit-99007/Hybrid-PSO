clear all;
close all;
clc
 %parameters
LB = [-10 -10];                          %parameter lower bound
UB = [10 10];                            %parameter upper bound
run = 10;
n= 100;                                %number of particle
Kmax = 50;                             %maximum iterations
r=rand(n,2);                           %random number
pos0=[];

pos0(:,1)= LB(1,1) + (UB(1,1)-LB(1,1))*r(:,1);
pos0(:,2)= LB(1,2) + (UB(1,2)-LB(1,2))*r(:,2);
p=pos0;                                %initial position
v = zeros(n,2);                        %initial velocity
k=1;
for r = 1:run
for k=1:Kmax
for i=1:n   %function output at position
    fv(i,1)=((p(i,1)^2+p(i,2)^2)/50-cos(p(i,1)/(i)^(1/2))*cos(p(i,2)/(i)^(1/2))+1 );
    %Updating pbest
    if k==1
        pbestval(i,1) = fv(i,1);
        pbest(i,:) = p(i,:);               %finding pbest
    elseif fv(i,1)< pbestval(i,1)          %checking condition if function value is less than pbest value
            pbestval(i,1)=fv(i,1);         %assign new pbest value as function value 
            pbest(i,:) = p(i,:);           %and correspondingly pbest will be at that position 
    end   
end

if k==1
    [gbestval(k,1),index] = min(fv);      %minimum f out and position 
    gbest(k,:) = pbest(index,:);          %finding gbest
elseif min(fv)< gbestval(k-1,1)           %if previous gbest value is greater than min of fv 
    [gbestval(k,1),index] = min(fv);      %updating gbestval to min of fv    
    gbest(k,:) = pbest(index,:);          %updating gbest to pbest
elseif min(fv) >= gbestval(k-1,1)         %if previous gbest value is less than min of fv
    gbest(k,:)=gbest(k-1,:);              %update gbest as previous one
    gbestval(k,1)=gbestval(k-1,1);        %update gbest value as previous gbest value
end
    
w=0.9 - (k/Kmax).*(0.9-0.4);              %calculation of w
c1=2.5 + (k/Kmax)*(0.5-2.5);              %calculation of c1
c2=0.5 + (k/Kmax)*(2.5-0.5);              %calculation of c2

%updating velocity and position
for i = 1:n
        v(i,1) = w*v(i,1) + c1*rand(1)*(pbest(i,1)-p(i,1)) + c2*rand(1)*(gbest(k,1)- p(i,1)); %updating velocity 1
        v(i,2) = w*v(i,1) + c1*rand(1)*(pbest(i,2)-p(i,2)) + c2*rand(1)*(gbest(k,2)- p(i,2)); %updating velocity 2
        p(i,:) = v(i,:) + p(i,:);         %updating positions
end
end
end
fprintf('PSO gbest val %10.20f \n',gbestval(k,1));
fprintf('PSO gbest %20.20f \n',gbest(k,1));

x=gbest(k,:);                           %Assuming x as gbest obtained
eps=0.000000000000001;

for k=1:Kmax
    %Defining Gradient Vector
    g(1,1) = (x(1)/25 + cos(x(2)/(2)^(1/2))*sin(x(1)));
    g(2,1) = (x(2)/25 + (cos((x(1))/(1)^(1/2)*sin((x(2))/(2)^(1/2))))/(1.414));
    %Defining Hessian Matrix
    H(1,1) = (cos((2^(1/2)*x(2))/2)*cos(x(1)) + 1/25 );
    H(1,2) = (-(2^(1/2)*sin((2^(1/2)*x(2))/2)*sin(x(1)))/2);
    H(2,1) = (-(2^(1/2)*sin((2^(1/2)*x(2))/2)*sin(x(1)))/2);
    H(2,2) = ((cos((2^(1/2)*x(2))/2)*cos(x(1)))/2 + 1/25);
    h = -inv(H)*g;                       %Calculating inverse of H and finally h
    if abs(g)<eps                        %checking condition that g value should not be less than eps
        break
    end
    xf(1)= x(1)+h(1);                   %updating x
    xf(2)= x(2)+h(2);                   %updating x
end
fv_N=(((xf(1)).^2+(xf(2)).^2)/50-cos((xf(1))/(i)^(1/2))*cos((xf(2))/(i)^(1/2))+1 );    %newton function value

fprintf('Newton gbest val %20.20f \n',fv_N);
fprintf('Newton gbest %20.20f \n',xf);
