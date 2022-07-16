%% CAAM 615 Final Project
%% Author: Javeria Ahmed

% Constants
dx = 10^-4;
As = 400*pi*10^-8;
R_a = 0.3; 
g_cl =1/15; C_m = 1; 


N = [250 7 249 1];
%  N(1) number of compartments in each daughter strand
%  N(2) 2N-1 number of comparments in junction where is the number of
%       daughter strands
%  N(3) number of compartments in mother
%  N(4) number of compartments in soma, usually one.

Ntemplate = [4 1 1 1];

%  Ntemp(1) number of daughter strands
%  Ntemp(2)  1 junction
%  Ntemp(3)  1 mother
%  Ntemp(4)  1 Soma

num_strand = sum(Ntemplate);
branch_count = Ntemplate(1);


a = 1e-4*ones(num_strand,1); 
A3 = 2*pi*a(end-1)*dx; rho = A3/As;
lambda = sqrt(a./(2*R_a*g_cl));

H =lgmdHines(N,Ntemplate,lambda,a, rho);
S = -(1/(dx^2))*H;

% D Transformation Matrix Construction
Db =[];
for i=1:branch_count
    Db = [Db a(i)*ones(1,N(1))];   
end

kk=0;
for i=(branch_count+1):(num_strand-1)
    Db = [Db a(i)*ones(1,N(2+kk))];
    kk=kk+1;
end
Db = [Db a(end-1)/rho];

D =   diag(Db);
sD =  diag((Db).^(1/2));
sDi = diag((Db).^(-1/2));

% Transformation to Symmetric Matrix
Sv = sD*S*sDi;

% Eigenvalue and Eigenvectors of Sv
[Q,Z] = eig(Sv);
z = diag(Z);
[sz,si] = sort(z);
Qs = Q(:,si);

x1 = N(end-1):(N(end-1)+N(1));
x2 = N(end-1):(N(end-1)+N(1));
x3 = 0:N(end-1);

for k=1:9
    q = sDi*Qs(:,k+1); %eigenvectors of H

    q1 = fliplr(q(1:N(1))'); %first daughter
    q2 = fliplr(q(((branch_count-1)*N(1)+1):branch_count*N(1))'); %last strand
    q3 = fliplr(q(branch_count*N(1)+N(2)+1:end)'); %mother

    subplot(3,3,k)
    plot(x3,q3,'r')
    hold on 
    plot(x2,[q3(end) q2],'k')
    plot(x1,[q3(end) q1],'r--')
    hold off 
  
end
xlabel('x/ \mum')
ylabel('w /cm^{-12}')
