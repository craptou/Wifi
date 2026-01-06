%%%%%%%
% main_compare.m
% TGn-B + ZF : comparaison non codé vs codé (hard/soft)
%%%%%%%
clear; close all; clc;

P = params_ofdm();
Nsym = 10000;

% Canal fixe pour comparaisons
[h, H, Heff] = tgnb_channel(P, P.seedChannel);

fprintf("Non code...\n");
TEB_unc = simulate_ofdm_tgnb(P, Nsym, h, Heff, false, 'hard');

fprintf("Code conv (hard)...\n");
TEB_hard = simulate_ofdm_tgnb(P, Nsym, h, Heff, true, 'hard');

fprintf("Code conv (soft)...\n");
TEB_soft = simulate_ofdm_tgnb(P, Nsym, h, Heff, true, 'soft');

figure;
semilogy(P.EbN0dB, TEB_unc,  "o-"); grid on; hold on;
semilogy(P.EbN0dB, TEB_hard, "s-");
semilogy(P.EbN0dB, TEB_soft, "d-");
xlabel("E_b/N_0 (dB)");
ylabel("TEB");
legend("ZF non codé", "ZF codé - Viterbi hard", "ZF codé - Viterbi soft", "Location","southwest");
title(sprintf("OFDM TGn-B (seed=%d) + AWGN - QPSK - ZF : hard vs soft", P.seedChannel));

figure;
subplot(2,1,1); plot(abs(H)); grid on; title("|H[k]| (TGn-B)");
subplot(2,1,2); plot(angle(H)); grid on; title("Phase(H[k]) (TGn-B)");

