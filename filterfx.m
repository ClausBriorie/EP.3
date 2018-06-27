function y=filterfx(b,a,x,B, sd)
% y = filterfx(b, a, x, B, sd)
% b = coeficientes do numerador
% a = coeficientes do denominador
% x = sinal de entrada
% B = número de bits
% sd = 0 se precisão simples, sd = 1 se precisão dupla

Nit = length(x);
y = zeros(1,Nit);
b = b(:)/a(1);
a = a/a(1);
a = a(:);
x = x(:)';
na = length(a);
nb = length(b);
if nargin < 5
    sd = 1;
end
if sd == 0
    for n = 1:Nit
        y(n) = 0;
        for k=1:min(n-1,na)
            y(n) = y(n) - quantize2(y(n-k)*a(k),B);
        end
        for k=1:min(n-1,nb)
            y(n) = y(n) + quantize2(x(n-k+1)*b(k),B);
        end
        y(n) = quantize2(y(n),B);
    end
else
    for n = 1:Nit
        y(n) = quantize2(-y([n-1:-1:max(1,n-na+1)])*a(2:min(n,na),1)+x([n:-1:max(1,n-nb+1)])*b(1:min(n,nb)),B);
    end
end