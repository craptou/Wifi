%%%%%%%
% main_compare.m
% OFDM sur canal TGn-B + AWGN, égalisation ZF
% Comparaison : non codé vs codé convolutif (Viterbi hard)
%%%%%%%

clear; close all; clc;

% ----- Paramètres -----
P = params_ofdm();
nb_symb_ofdm = 10000;     % à réduire si c'est trop long
fprintf("Nfft=%d, Nu=%d, Ncp=%d, QPSK(M=%d), nbSymb=%d\n", ...
    P.Nfft, P.Nu, P.Ncp, P.M, nb_symb_ofdm);


% ----- Canal fixe (même réalisation pour comparer proprement) -----
[h, H, Heff] = tgnb_channel(P, P.seed_canal);

% ----- Simulations -----
fprintf("Simulation non codee...\n");
TEB_non_code = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, false);

fprintf("Simulation codee (conv + Viterbi hard)...\n");
TEB_code_hard = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, true);

% ----- Tracé TEB -----
figure;
semilogy(P.EbN0_dB, TEB_non_code, 'o-'); grid on; hold on;
semilogy(P.EbN0_dB, TEB_code_hard, 's-');
xlabel("E_b/N_0 (dB)");
ylabel("TEB");
legend("ZF non codé", "ZF codé (Viterbi hard)", "Location","southwest");
title(sprintf("OFDM TGn-B (seed=%d) + AWGN - QPSK - Egalisation ZF", P.seed_canal));

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
