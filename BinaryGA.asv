% Binary-coded GA
%---------------------------------------------------------
clc;
clear all;
close all;

pop=input('Enter population size: ');      %population size
str=10;                                    %string length
pc=input('Enter crossover probability: '); %crossover probability
pm=input('Enter mutation probability: ');  %mutation probability
x1_min=0;                                  %minimum value of x1
x1_max=0.5;                                %maximum value of x2
x2_min=0;                                  %minimum value of x1
x2_max=0.5;                                %maximum value of x2

itr=0;                                     %iteration
x=randi([0 1],pop,str);                    %creation of first generation
while (itr<=500)                            
for i=1:pop                                %dividing the string length into x1 and x2       
    for r=1:str
        if (r<=str/2)
            x1(i,r)=x(i,r);                %matrix of substring of x1 
        else
            x2(i,r-(str/2))=x(i,r);        %matrix of substring of x2
        end
    end
end

dvx1=zeros(pop,1);
dvx2=zeros(pop,1);

% Decoding
for i=1:pop
    for r=1:str/2
        dvx1(i,1)=dvx1(i,1)+x1(i,r)*2^(str/2-r);
        dvx2(i,1)=dvx2(i,1)+x2(i,r)*2^(str/2-r);
    end
end 

x1_val=zeros(pop,1);
x2_val=zeros(pop,1);

for i=1:pop
    x1_val(i,1)=x1_val(i,1)+x1_min+(((x1_max-x1_min)/(2^(str/2)-1))*dvx1(i,1));
    x2_val(i,1)=x2_val(i,1)+x2_min+(((x2_max-x2_min)/(2^(str/2)-1))*dvx2(i,1));
end


x1_real(itr+1)=mean(x1_val(:,1));
x2_real(itr+1)=mean(x2_val(:,1));


% Fitness calculation
fitx=zeros(pop,1);
for i=1:pop
    fitx(i,1)=func(x1_val(i,1),x2_val(i,1));
end

min_fitness(itr+1)=min(fitx);
max_fitness(itr+1)=max(fitx);
avg_fitness(itr+1)=sum(fitx)/pop;

fitnew=fitx/sum(fitx);

% Selection (Roulette Wheel Selection)
for i=2:pop
    fitnew(i)=fitnew(i)+fitnew(i-1);
end

for i=1:pop                                     % Mating pool generation
    random=rand;
    for r=1:pop
        if random<=fitnew(r)
            matingpool(i,:)=x(r,:);
            break
        end
    end
end
    
    randomx=randperm(pop,pop);
    for i=1:length(randomx)
        newmat(i,:)= matingpool(randomx(i),:);
    end
        
   
% Single-point crossover    
    for i=1:2:pop
        randomr=rand;
    if randomr<pc
               y1 = newmat(i,:);
               y2 = newmat(i+1,:);
               r = randi([1, str-1]);
    
               m1 = [y1(1:r) y2(r+1:str)];
               m2 = [y2(1:r) y1(r+1:str)];
               
               newmat(i,:) = m1;
               newmat(i+1,:) = m2;
    end
    end

% Mutation
for i=1:pop
    for r=1:str
        randomz=rand;
        if randomz<pm
            if newmat(i,r)==0
               newmat(i,r)=1;
            else
               newmat(i,r) = 0;
            end
        end
    end
end

x=newmat;
itr=itr+1;

end

%% Plotting
figure(1)
plot(1:itr,avg_fitness);
axis([1 itr 0.5 1]);
xlabel('Generations');
ylabel('Average fitness');
legend('Average fitness');
title('Average fitness vs No. of generations');
hold off;


figure(2)
plot(1:itr,max_fitness);
axis([1 itr 0.5 1]);
xlabel('Generations');
ylabel('Fitness');
hold on;

plot(1:itr,min_fitness);
axis([1 itr 0.5 1]);
xlabel('Generations');
ylabel('Fitness');
legend('Maximum fitness','Minimum fitness');
title('Minimum fitness and Maximum fitness vs No. of generations');
hold off;

figure(3)
plot(1:itr,x1_real);
axis([1 itr 0 1]);
xlabel('Generations');
ylabel('Variable values');
hold on;

plot(1:itr,x2_real);
axis([1 itr 0 1]);
xlabel('Generations');
ylabel('Variable values');
legend('x1','x2');
title('Optimal Solution');
hold off;

minfunc_value=(1/max(fitx)-1)^0.5
x1_value=min(x1_val)
x2_value=min(x2_val)


%% Objective function
function f=func(x1,x2)

eqn=fileread('Input.txt');
fh=str2func(eqn);
f1=fh(x1,x2);
f=1/(1+f1^2);
end