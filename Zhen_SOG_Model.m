%Stochastic Optimal Growth with Capital Persistence. Homework for Zhen.
clear

global v0 beta delta alpha k_m k0 M a0 j;

alpha = 0.3;
beta = 0.96;
delta = 0.1; % depreciation rate (annual)
epsilon = 0.05;
z_H = 1 + epsilon;
z_L = 1 - epsilon;
amat = [z_H z_L]';
p = 0.5;
M = zeros(2,2); %M for markov transition kernel
M(1,1) = p;
M(2,2) = p;
M(2,1) = 1-p;
M(1,2) = 1-p;

tol = 0.01;
maxits = 1000;
dif = tol+1000;
its = 0;

kgrid = 49; % grid points + 1


kmin = 0.5;
kmax = 7;
grid = (kmax-kmin)/kgrid;

k_m = kmin:grid:kmax;
k_m = k_m';
[N,n] = size(k_m);

polfun = zeros(kgrid+1,2);

v0 = zeros(N,2);
dif = 10;
its = 0;
k11 = zeros(N,2);
while dif > tol %& its < maxits
    
    for j = 1:2
        for i = 1:N
            k0 = k_m(i,1);
            a0 = amat(j,1);
            k1 = fminbnd(@valfun_stoch,kmin,kmax);
            v1(i,j) = -valfun_stoch(k1);
            k11(i,j) = k1;
        end
    end
    
    dif = norm(v1-v0);
    v0 = v1;
    its = its+1
end

for i = 1:N
    con(i,1) = k_m(i,1)^(alpha) - k11(i,1) + (1-delta)*k_m(i,1);
    polfun(i,:) = k_m(i,:)^(alpha) - k11(i,:) + (1-delta)*k_m(i,:);
end



figure
plot(k_m,polfun,'LineWidth', 1)
xlabel('k')
ylabel('c')

figure
plot(k_m, [k11, k_m], 'Linewidth',1)
xlabel('k')
ylabel('k next')
