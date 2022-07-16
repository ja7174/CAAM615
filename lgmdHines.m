%% CAAM 615 Final Project
%% Author: Javeria Ahmed

function H = lgmdHines(N,Ntemplate,lambda,a, rho)

total_len = N*Ntemplate';
H = zeros(total_len,total_len);
branch_count = Ntemplate(1);
branch_len = N(1);
r = a./a(branch_count+1);

 %% Daughter Branches
r1 = 1; r2 = 0;
for i=1:branch_count
   r2 = r2+ branch_len;%N(1,i);
   H_sub = HinesSub(branch_len);%
   H(r1:r2,r1:r2) = (lambda(i)^2)* H_sub;%
   r1 = r2 +1;
end
junction_index = N(1)*Ntemplate(1);
for i=1:branch_count
    H(i*branch_len,2*i-1+junction_index) = (lambda(i)^2);  %%*       
end
 
 %% Junction 

junction_len = 2*Ntemplate(1)-1;
junction_sub = zeros(junction_len,junction_len);
m = N(1)*Ntemplate(1)+junction_len+1;

for i=1:junction_len %this iterates over the rows
    kk = branch_count+1;
    if mod(i,2)==1 % bridge
        if i==1 %start
            junction_sub(i,1) = -lambda(kk)^2; 
            junction_sub(i,2) = lambda(kk)^2;
        elseif i==junction_len %end
            junction_sub(i,junction_len) = -1*lambda(kk)^2;
            junction_sub(i,junction_len-1) = 1*lambda(kk)^2;            
        else %middle
            junction_sub(i,i+1) = 1*lambda(kk)^2;
            junction_sub(i,i-1) = 1*lambda(kk)^2;
            junction_sub(i,i) = -2*lambda(kk)^2;
        end
        j = (i+1)/2;
        junction_sub(i,i) = junction_sub(i,i) -r(j)*lambda(j)^2;
        

    elseif i==Ntemplate(1) % central junction 
            
            junction_sub(i,i+1) = 1*lambda(kk)^2;
            junction_sub(i,i-1) = 1*lambda(kk)^2;
            junction_sub(i,i) = -2*lambda(kk)^2;
            junction_sub(i,i) = junction_sub(i,i)-r(end)*lambda(end)^2;
            H(junction_index+i,m)= -r(end)*lambda(end)^2;
            
            
    else %connector
            
            junction_sub(i,i+1) = 1*lambda(kk)^2;
            junction_sub(i,i-1) = 1*lambda(kk)^2;
            junction_sub(i,i) = -2*lambda(kk)^2;

    end
end

for i=1:branch_count
    H(2*i-1+junction_index,i*branch_len) = r(i)*(lambda(i)^2);  %%       
end

rj_start = junction_index+1;
rj_end= junction_index+1+junction_len -1;
H(rj_start:rj_end,rj_start:rj_end) = junction_sub;

 %% Mother
% 
mother_len = N(end-1);
M_sub = HinesSub(mother_len);
M_sub(1,1) =-2;
M_sub(end,end) =-2;

rm_start = m;
rm_end = m+mother_len-1;

H(rm_start:rm_end,rm_start:rm_end) = (lambda(end-1)^2)*M_sub; 
H(m,junction_index+Ntemplate(1))= -lambda(end-1)^2; 
H(end-1,end) = lambda(end-1)^2;
 %% Soma
 H(end,end-1)=rho*lambda(end)^2;
 H(end,end)=-rho*lambda(end)^2;

end