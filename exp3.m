% Parâmetros
f0 = 300;
fa = 40e3;
Ta = 1/fa;
dur = 3;

% Vetor com os tempos amostrados
vec_t_am = 0:Ta:dur;

% Sinal
s = 0.45*cos(2 * pi * f0 * vec_t_am) + 0.45*cos(6 * pi * f0 * vec_t_am);

[d, tamanho_s] = size(s);

% Ouvindo sinal criado
% sound(s, fa);

% Criando sinal quantizado com 5 bits
B = 5;
sq = quantize2(s, B);

%% Ouvindo sinal com quantização
% é possí­vel perceber o ruído introduzido pela quantização
% sound(s, fa);

%% Calculando relação sinal-ruído:
% Potência do sinal original

pot_s = 0;
for n = 1:tamanho_s
    aux = s(n)^2;
    pot_s = pot_s + aux;
end
pot_s = pot_s/tamanho_s;
pot_s_dB = -20*log10(pot_s);

ruido = s - sq;
% Potência do ruído
pot_ruido = 0;
for n = 1:tamanho_s
    aux = ruido(n)^2;
    pot_ruido = pot_ruido + aux;
end
pot_ruido = pot_ruido/tamanho_s;
pot_ruido_dB = 20*log10(pot_ruido);

SNR = 20*log10(pot_s/pot_ruido);

%% Projetando um filtro FIR
%       Especificações do filtro:

% Faixa de passagem [rad/amostra] :
f_pass = 0.05*pi*fa/(2*pi);
% Faixa de rejeição [rad/amostra] :
f_rej  = 0.2*pi*fa/(2*pi);

% Vetor de frequências
f = [f_pass f_rej];

% Atenuação máxima na faixa de rejeição ([dB] e linear):
at_minRej_dB = 40;
at_minRej = 10^(-at_minRej_dB/20);

% Queda máxima na faixa de passagem:
q_max = 0.05; % = delta_p
at_maxPas_dB = -20*log10(1-q_max);
at_maxPas = 10^(-at_maxPas_dB/20);

% Vetor de ganhos (indicando que é um filtro passa-baixas)
a = [1 0];
dev = [q_max, at_minRej];

% Projeto dos filtros
cell_FIR = firpmord(f, a, dev, fa, 'cell');

% Ajuste no número de coeficientes (pois as especificações não foram
% atingidas usando os valores estimados por firpmord)
cell_FIR{1} = cell_FIR{1}+6;

h_F = firpm(cell_FIR{:});
figure(1)
freqz(h_F, 1);
title(sprintf('Filtro Parks-McClellan passa-baixas (n = %d)', cell_FIR{1}))

%% Projetando um filtro IIR
Wp = f_pass/(fa/2);
Ws = f_rej/(fa/2);

[n, Wp_2] = ellipord(Wp, Ws, at_maxPas_dB, at_minRej_dB);

[z_iir, p_iir] = ellip(n, at_maxPas_dB, at_minRej_dB, Wp_2);

h_I = impz(z_iir,p_iir);

figure(2)
freqz(z_iir, p_iir);
title(sprintf('Filtro ElÃ­ptico passa-baixas (n = %d)', n))

%% Exercício 1: função valor esperado de cada filtro

% Somatórios dos coeficientes
somatorio_FIR = sum(h_F);
somatorio_IIR = sum(h_I);
somatorio_sq_FIR = sum(h_F.^2);
somatorio_sq_IIR = sum(h_I.^2);

%% Exercício 2:
Y_F = filter(h_F, 1, sq);        % saída do FIR para entrada quantizada
Y_I = filter(z_iir, p_iir, sq);  % saída do IIR para entrada quantizada

y_F = filter(h_F, 1, s);        % saída do FIR para entrada pura
y_I = filter(z_iir, p_iir, s);  % saída do IIR para entrada pura

e_F = Y_F - y_F;
e_I = Y_I - y_I;

% Cálculos dos valores DC
valor_DC_FIR = mean(e_F);
valor_DC_IIF = mean(e_I);

% Cálculos das potências médias
PotMedia_teorica_ruido_FIR = ((2^(-2*B))/3)*somatorio_sq_FIR;
PotMedia_teorica_ruido_IIR = ((2^(-2*B))/3)*somatorio_sq_IIR;

PotMedia_exp_sinalQuantizado_FIR = 0;
for n = 1:tamanho_s
    aux = Y_F(n)^2;
    PotMedia_exp_sinalQuantizado_FIR = PotMedia_exp_sinalQuantizado_FIR + aux;
end
PotMedia_exp_sinalQuantizado_FIR = PotMedia_exp_sinalQuantizado_FIR/tamanho_s

PotMedia_exp_sinalQuantizado_IIR = 0;
for n = 1:tamanho_s
    aux = Y_I(n)^2;
    PotMedia_exp_sinalQuantizado_IIR = PotMedia_exp_sinalQuantizado_IIR + aux;
end
PotMedia_exp_sinalQuantizado_IIR = PotMedia_exp_sinalQuantizado_IIR/tamanho_s

% Cálculos de SNRs nas saídas
SNR_teorico_FIR = PotMedia_exp_sinalQuantizado_FIR/PotMedia_teorica_ruido_FIR
SNR_teorico_IIR = PotMedia_exp_sinalQuantizado_IIR/PotMedia_teorica_ruido_IIR

%% Exercício 3
% Diferença entre saídas para entradas pura e quantizada
y_F = filter(h_F, 1, s);        % saída do FIR para entrada pura
y_I = filter(z_iir, p_iir, s);  % saída do IIR para entrada pura

e_F = Y_F - y_F;
e_I = Y_I - y_I;

PotMedia_exp_ruido_FIR = 0;
for n = 1:tamanho_s
    aux = e_F(n)^2;
    PotMedia_exp_ruido_FIR = PotMedia_exp_ruido_FIR + aux;
end
PotMedia_exp_ruido_FIR = PotMedia_exp_ruido_FIR/tamanho_s

PotMedia_exp_ruido_IIR = 0;
for n = 1:tamanho_s
    aux = e_I(n)^2;
    PotMedia_exp_ruido_IIR = PotMedia_exp_ruido_IIR + aux;
end
PotMedia_exp_ruido_IIR = PotMedia_exp_ruido_IIR/tamanho_s

% É possível perceber que os valores das potências dos ruídos são próximos

%% Exercício 4
Y_F2 = filterfx(h_F, 1, sq, 2*B, 1);        % saída do FIR para entrada quantizada, considerando a introdução de erros de quantização
Y_I2 = filterfx(z_iir, p_iir, sq, 2*B, 1);  % saída do IIR para entrada quantizada, considerando a introdução de erros de quantização

%% Exercício 5
e_F2 = Y_F2 - y_F;
e_I2 = Y_I2 - y_I;

PotMedia_exp_fx_ruido_FIR = 0;
for n = 1:tamanho_s
    aux = e_F2(n)^2;
    PotMedia_exp_fx_ruido_FIR = PotMedia_exp_fx_ruido_FIR + aux;
end
PotMedia_exp_fx_ruido_FIR = PotMedia_exp_fx_ruido_FIR/tamanho_s

PotMedia_exp_fx_ruido_IIR = 0;
for n = 1:tamanho_s
    aux = e_I2(n)^2;
    PotMedia_exp_fx_ruido_IIR = PotMedia_exp_fx_ruido_IIR + aux;
end
PotMedia_exp_fx_ruido_IIR = PotMedia_exp_fx_ruido_IIR/tamanho_s

SNR_fx_FIR = PotMedia_exp_sinalQuantizado_FIR/PotMedia_exp_fx_ruido_FIR
SNR_fx_IIR = PotMedia_exp_sinalQuantizado_IIR/PotMedia_exp_fx_ruido_IIR
