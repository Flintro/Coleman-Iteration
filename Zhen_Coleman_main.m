clear
global beta alpha sigma k_grid k0 delta policy_guess M amat a0 j;
alpha = 0.3;
beta = 0.96;
sigma = 0.9;
delta = 0.1;

%Stochasticity parameters
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


N = 100;
kmin = 0.5;
kmax = 7;
k_grid = linspace(kmin,kmax,N);
policy_guess = zeros(2,N);
policy_guess(1,:) = 0.3*k_grid.^alpha;
policy_guess(2,:) = 0.3*k_grid.^alpha;
%policy_guess = zeros(1,N);
k_new = zeros(2,N);

maxits = 20;
its = 0;

while its < maxits
    for j = 1:2
        for i = 1:N
            k0 = k_grid(1,i);
            a0 = amat(j,1);
            %k1 = fminbnd(@qe_euler_diff, kmin, kmax);
            k1 = fsolve( @(k) euler_diff(k) , 1 );
            k_new(j,i) = k1;
        end
        
    
    %policy_guess(j,:) = k_grid.^alpha + (1- delta)*k_grid - k_new(j,:);
    
    end
    policy_guess = a0*k_grid.^alpha + (1- delta)*k_grid - k_new;
    its = its+1
    
end

figure
plot(k_grid, k_new, 'Linewidth', 1.5)
xlabel('k init');
ylabel('k next')
hold on
plot(k_grid, k_grid)
hold off