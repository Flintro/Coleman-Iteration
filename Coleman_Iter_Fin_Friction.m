%Solve the Household problem with borrowing constraint
%using Coleman Policy iteration
clear
global beta r bbar a y a_grid policy_guess M j;
beta = 0.96; %discount
bbar = 0; %-0.01; %borrowing constraint
r = 0.02; %interest rate

%Stochasticity parameters
epsilon = 0.15;
y_H = 1 + epsilon;
y_L = 1 - epsilon;
amat = [y_H y_L]';
p = 0.8;
M = zeros(2,2); %M for markov transition kernel
M(1,1) = p;
M(2,2) = p;
M(2,1) = 1-p;
M(1,2) = 1-p;


N = 100;%grid size
amin = 0; %grid bounds
amax = 3;
a_grid = linspace(amin,amax,N);
policy_guess = zeros(2,N);
policy_guess(1,:) = 0.5*a_grid;
policy_guess(2,:) = 0.5*a_grid;
policy_new = repmat(policy_guess,1);
a_new = zeros(2,N);

tol = 1e-06;
maxits = 20;
its = 0;

while its < maxits %& max(abs((policy_guess - policy_guess_old))) > tol
    for j = 1:2
        for i = 1:N
            a = a_grid(1,i);
            y = amat(j,1);
            %k1 = fminbnd(@abs(euler_diff_test), kmin, kmax);
            %a1 = fsolve( @(a1) euler_diff_test(a1), 0.5 ); %try 0.5, 1 and 2 for guess point
            a1 = (1+r)*a + y - policy_guess(j,i);
            if a1>= bbar
                c1 = 0;
                for q = 1:2
                    policy_func = interp1(a_grid', policy_guess(q,:), a1, 'linear');
                    c1 = c1 + beta*(1+r)*(1 ./ policy_func)*M(j,q)';        %1/c_t+1
                end
                policy_new(j,i) = 1/c1;
            else
                policy_new(j,i) = (1+r)*a_grid(:,i) + y - bbar;
            end
            
            a_new(j,i) = max(a1, bbar);
            
        end
        
    
    
    end
    %Update the consumption policy, 
    %BUT, note we use partial weighting to assure convergence stability
    policy_guess = 0.8 .* policy_new + 0.2 .* policy_guess;
    
    its = its+1
    
end

figure
plot(a_grid, a_new, 'Linewidth', 1)
xlabel('A init');
ylabel('A next')
hold on
plot(a_grid, a_grid)
hold off

% figure
% plot(a_grid, policy_guess, 'Linewidth', 1.5)
% xlabel('A init');
% ylabel('Optimal Consumption')
% hold on
% plot(a_grid, a_grid)
% hold off