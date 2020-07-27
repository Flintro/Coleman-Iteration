function val=valfun_stoch(k)

% This program gets the value function for a stochastic growth model with
% CRRA utility

global v0 beta delta alpha k_m k0 M a0 s j

g = interp1(k_m,v0,k,'linear');

c = a0*k0^alpha - k + (1-delta)*k0; % consumption
if c<=0
    val = -8888888888888888-800*abs(c); % keeps it from going negative
else

val = log(c) + beta*(g*M(j,:)');
end
val = -val;