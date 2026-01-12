% main_entrelacement.m - Impact de l'entrelacement

clear; close all; clc;

P = params_ofdm();
nb_symb_ofdm = 1000;

[h, H, Heff] = tgnb_channel(P, P.seed_canal);

fprintf("Non code...\n");
TEB_non_code = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'none', 0);

fprintf("R=1/2 sans entrelacement...\n");
P.activer_entrelacement = false;
TEB_code_1_2 = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 1/2);

fprintf("R=1/2 avec entrelacement...\n");
P.activer_entrelacement = true;
TEB_code_1_2_int = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 1/2);

fprintf("R=2/3 sans entrelacement...\n");
P.activer_entrelacement = false;
TEB_code_2_3 = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 2/3);

fprintf("R=2/3 avec entrelacement...\n");
P.activer_entrelacement = true;
TEB_code_2_3_int = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 2/3);

fprintf("R=3/4 sans entrelacement...\n");
P.activer_entrelacement = false;
TEB_code_3_4 = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 3/4);

fprintf("R=3/4 avec entrelacement...\n");
P.activer_entrelacement = true;
TEB_code_3_4_int = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, 'conv', 3/4);

% Tracé
figure;

subplot(1,3,1);
semilogy(P.EbN0_dB, TEB_non_code, 'k:o'); grid on; hold on;
semilogy(P.EbN0_dB, TEB_code_1_2, 'b-s');
semilogy(P.EbN0_dB, TEB_code_1_2_int, 'r-d');
xlabel("E_b/N_0 (dB)"); ylabel("TEB");
legend("Non codé", "R=1/2 sans", "R=1/2 avec", "Location", "southwest");
title("R=1/2");

subplot(1,3,2);
semilogy(P.EbN0_dB, TEB_non_code, 'k:o'); grid on; hold on;
semilogy(P.EbN0_dB, TEB_code_2_3, 'b-s');
semilogy(P.EbN0_dB, TEB_code_2_3_int, 'r-d');
xlabel("E_b/N_0 (dB)"); ylabel("TEB");
legend("Non codé", "R=2/3 sans", "R=2/3 avec", "Location", "southwest");
title("R=2/3");

subplot(1,3,3);
semilogy(P.EbN0_dB, TEB_non_code, 'k:o'); grid on; hold on;
semilogy(P.EbN0_dB, TEB_code_3_4, 'b-s');
semilogy(P.EbN0_dB, TEB_code_3_4_int, 'r-d');
xlabel("E_b/N_0 (dB)"); ylabel("TEB");
legend("Non codé", "R=3/4 sans", "R=3/4 avec", "Location", "southwest");
title("R=3/4");
