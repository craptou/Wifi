% main_compare.m - Comparaison non code vs Conv vs LDPC

clear; close all; clc;

P = params_ofdm();
nb_symb_ofdm = 1000;

[h, H, Heff] = tgnb_channel(P, P.seed_canal);

fprintf("Simulation non codee...\n");
TEB_non_code = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'none', 0);

fprintf("Simulation Conv R=1/2...\n");
P.activer_entrelacement = true;
TEB_conv_1_2 = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 1/2);

fprintf("Simulation LDPC R=1/2...\n");
TEB_ldpc = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'ldpc', 0);

% Tracé
figure;
semilogy(P.EbN0_dB, TEB_non_code, 'k:o', 'LineWidth', 2); 
grid on; hold on;
semilogy(P.EbN0_dB, TEB_conv_1_2, 'b-s', 'LineWidth', 2);
semilogy(P.EbN0_dB, TEB_ldpc, 'r-d', 'LineWidth', 2);
xlabel("E_b/N_0 (dB)");
ylabel("TEB");
legend("Non codé", "Conv R=1/2", "LDPC R=1/2", "Location", "southwest");
title("Comparaison Codage Conv vs LDPC");

% Canal
figure;
subplot(2,1,1);
plot(abs(H)); grid on;
title("|H[k]|");
xlabel("Indice k");

subplot(2,1,2);
plot(angle(H)); grid on;
title("Phase(H[k])");
xlabel("Indice k");
