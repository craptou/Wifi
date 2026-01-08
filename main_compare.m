%%%%%%%
% main_compare.m
% OFDM sur canal TGn-B + AWGN, égalisation ZF
% Comparaison : non codé vs codé convolutif vs codé LDPC
%%%%%%%

clear; close all; clc;

% ----- Paramètres -----
P = params_ofdm();
nb_symb_ofdm = 1000;     % réduit pour tests rapides (augmenter pour résultats finaux)
fprintf("Nfft=%d, Nu=%d, Ncp=%d, QPSK(M=%d), nbSymb=%d\n", ...
    P.Nfft, P.Nu, P.Ncp, P.M, nb_symb_ofdm);

% ----- Canal fixe (même réalisation pour comparer proprement) -----
[h, H, Heff] = tgnb_channel(P, P.seed_canal);

% ----- Simulations : comparaison LDPC vs Convolutif -----
fprintf("\n===== Comparaison Non code vs Conv vs LDPC =====\n");

fprintf("Simulation non codee...\n");
TEB_non_code = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'none', 0);

fprintf("Simulation codee convolutif R=1/2...\n");
P.activer_entrelacement = false;
TEB_conv_1_2 = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 1/2);

fprintf("Simulation codee LDPC R=1/2...\n");
TEB_ldpc = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'ldpc', 0);

% ----- Tracé TEB : Comparaison Conv vs LDPC -----
figure('Position', [100 100 1000 700]);

semilogy(P.EbN0_dB, TEB_non_code, 'k:o', 'LineWidth', 2, 'MarkerSize', 8); 
grid on; hold on;
semilogy(P.EbN0_dB, TEB_conv_1_2, 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(P.EbN0_dB, TEB_ldpc, 'r-d', 'LineWidth', 2, 'MarkerSize', 8);

xlabel("E_b/N_0 (dB)", 'FontSize', 12);
ylabel("TEB", 'FontSize', 12);
legend("Non codé", "Conv R=1/2", "LDPC R=1/2", "Location","southwest", 'FontSize', 11);
title(sprintf("OFDM TGn-B - Comparaison Codage Conv vs LDPC (seed=%d)", P.seed_canal), ...
    'FontSize', 13);
grid on;

% Ajout d'annotations
text(5, 1e-1, sprintf('LDPC: N=%d, K=%d, iter=%d', P.ldpc_N, P.ldpc_K, P.ldpc_nb_iter), ...
    'FontSize', 10, 'BackgroundColor', 'white');

fprintf("\n===== Simulations terminees =====\n");

% ----- Visualisation du canal -----
figure;
subplot(2,1,1);
plot(abs(H)); grid on;
title("|H[k]| (canal TGn-B)");
xlabel("Indice sous-porteuse k");
ylabel("|H[k]|");

subplot(2,1,2);
plot(angle(H)); grid on;
title("Phase(H[k]) (canal TGn-B)");
xlabel("Indice sous-porteuse k");
ylabel("angle(H[k])");
