
function y = quantize2(x,B)  
% y = quantize2(x,B) 
% The function rounds x into a binary fixed point 
% representation with B bits, including the sign bit. 

% If overflow occurs, the largest possible 
% fixed-point representation is returned. 
 
% The largest numbers that can be represented are
% + 2^(k-1)-1 or - 2^(k-1)

% If x is less than these extreme values, the routine finds
% the number of bits that are needed to represent the integer
% part of x and uses the remaining bits for its fractional part.

y = 0;
B1 = B-1;
fac = 2^(B1);
invfac = 2^(-B1);
xint = round(x*fac);
xauxpos = xint >= (fac - 1); % positive overflow
xauxneg = xint <= -fac; % negative overflow
y = invfac * (xint .* (~xauxpos) .* (~xauxneg))+ (1-invfac)*xauxpos - xauxneg; % com overflow
%y = invfac * xint; % sem overflow

%y=x;