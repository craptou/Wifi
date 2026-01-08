%%%%%%%
% main_entrelacement.m
% OFDM sur canal TGn-B + AWGN, égalisation ZF
% Comparaison : impact de l'entrelacement sur différents rendements
%%%%%%%

clear; close all; clc;

% ----- Paramètres -----
P = params_ofdm();
nb_symb_ofdm = 1000;     % réduit pour tests rapides
fprintf("Nfft=%d, Nu=%d, Ncp=%d, QPSK(M=%d), nbSymb=%d\n", ...
    P.Nfft, P.Nu, P.Ncp, P.M, nb_symb_ofdm);

% ----- Canal fixe (même réalisation pour comparer proprement) -----
[h, H, Heff] = tgnb_channel(P, P.seed_canal);

% ----- Simulations : impact entrelacement -----
fprintf("\n===== Comparaison avec/sans entrelacement =====\n");

fprintf("Simulation non codee...\n");
TEB_non_code = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'none', 0);

fprintf("Simulation R=1/2 SANS entrelacement...\n");
P.activer_entrelacement = false;
TEB_code_1_2 = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 1/2);

fprintf("Simulation R=1/2 AVEC entrelacement...\n");
P.activer_entrelacement = true;
TEB_code_1_2_int = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 1/2);

fprintf("Simulation R=2/3 SANS entrelacement...\n");
P.activer_entrelacement = false;
TEB_code_2_3 = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 2/3);

fprintf("Simulation R=2/3 AVEC entrelacement...\n");
P.activer_entrelacement = true;
TEB_code_2_3_int = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 2/3);

fprintf("Simulation R=3/4 SANS entrelacement...\n");
P.activer_entrelacement = false;
TEB_code_3_4 = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 3/4);

fprintf("Simulation R=3/4 AVEC entrelacement...\n");
P.activer_entrelacement = true;
TEB_code_3_4_int = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 3/4);

% ----- Tracé TEB : Comparaison avec/sans entrelacement -----
figure('Position', [100 100 1400 600]);

% Subplot 1 : R=1/2
subplot(1,3,1);
semilogy(P.EbN0_dB, TEB_non_code, 'k:o', 'LineWidth', 1.5); grid on; hold on;
semilogy(P.EbN0_dB, TEB_code_1_2, 'b-s', 'LineWidth', 1.5);
semilogy(P.EbN0_dB, TEB_code_1_2_int, 'r-d', 'LineWidth', 1.5);
xlabel("E_b/N_0 (dB)", 'FontSize', 11); ylabel("TEB", 'FontSize', 11);
legend("Non codé", "R=1/2 sans int.", "R=1/2 avec int.", "Location","southwest");
title("Code Conv R=1/2 - Impact entrelacement", 'FontSize', 12);

% Subplot 2 : R=2/3
subplot(1,3,2);
semilogy(P.EbN0_dB, TEB_non_code, 'k:o', 'LineWidth', 1.5); grid on; hold on;
semilogy(P.EbN0_dB, TEB_code_2_3, 'b-s', 'LineWidth', 1.5);
semilogy(P.EbN0_dB, TEB_code_2_3_int, 'r-d', 'LineWidth', 1.5);
xlabel("E_b/N_0 (dB)", 'FontSize', 11); ylabel("TEB", 'FontSize', 11);
legend("Non codé", "R=2/3 sans int.", "R=2/3 avec int.", "Location","southwest");
title("Code Conv R=2/3 - Impact entrelacement", 'FontSize', 12);

% Subplot 3 : R=3/4
subplot(1,3,3);
semilogy(P.EbN0_dB, TEB_non_code, 'k:o', 'LineWidth', 1.5); grid on; hold on;
semilogy(P.EbN0_dB, TEB_code_3_4, 'b-s', 'LineWidth', 1.5);
semilogy(P.EbN0_dB, TEB_code_3_4_int, 'r-d', 'LineWidth', 1.5);
xlabel("E_b/N_0 (dB)", 'FontSize', 11); ylabel("TEB", 'FontSize', 11);
legend("Non codé", "R=3/4 sans int.", "R=3/4 avec int.", "Location","southwest");
title("Code Conv R=3/4 - Impact entrelacement", 'FontSize', 12);

sgtitle(sprintf("OFDM TGn-B - Entrelacement IEEE 802.11a (seed=%d)", P.seed_canal), ...
    'FontSize', 13);

fprintf("\n===== Simulations terminees =====\n");
