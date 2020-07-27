function diff = euler_diff(k)

global beta alpha sigma k_grid k0 delta policy_guess M amat a0 j;

c = a0*k0^alpha + (1-delta)*k0 - k;
if c<=1e-02
    
    diff = 888888888888888888888;
%     c = 0;
%     k = k0^alpha + (1-delta)*k0;
else
policy_func = interp1(k_grid', policy_guess', k);%(j,:), k);

diff = ( 1/c - beta*(1 ./ policy_func)*M(j,:)'*(a0*alpha*k^(alpha-1) + (1 - delta) ) ); %still need to multiply by markov kernel
end



end