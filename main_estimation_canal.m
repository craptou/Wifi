%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_estimation_canal.m
% PARTIE 8 : Estimation du canal de propagation par porteuses pilotes
%
% Compare les performances TEB entre :
%   - Canal parfaitement connu (référence)
%   - Canal estimé à partir des 4 porteuses pilotes IEEE 802.11a
%
% Méthode d'estimation :
%   1. Estimation LS sur les 4 pilotes : Hp[k] = Yp[k] / Xp[k]
%   2. Interpolation linéaire vers les sous-porteuses de données
%   3. Moyennage temporel (canal stationnaire)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

% ----- Paramètres -----
P = params_ofdm();
nb_symb_ofdm = 1000;  % nombre de symboles OFDM

fprintf("===== PARTIE 8 : Estimation du canal par pilotes =====\n");
fprintf("Nfft=%d, Nu=%d, Ncp=%d, QPSK(M=%d), nbSymb=%d\n", ...
    P.Nfft, P.Nu, P.Ncp, P.M, nb_symb_ofdm);
fprintf("Pilotes aux indices : [%s]\n", num2str(P.idx_pilotes));
fprintf("Symboles pilotes : [%s]\n", num2str(P.symboles_pilotes'));

% ----- Canal fixe (même réalisation pour comparaison équitable) -----
[h, H, Heff] = tgnb_channel(P, P.seed_canal);

% ----- Simulation avec canal parfait -----
fprintf("\nSimulation avec canal PARFAIT...\n");
TEB_parfait = simulate_ofdm_pilotes(P, nb_symb_ofdm, h, Heff, 'parfait');

% ----- Simulation avec estimation par pilotes -----
fprintf("Simulation avec estimation par PILOTES...\n");
TEB_pilotes = simulate_ofdm_pilotes(P, nb_symb_ofdm, h, Heff, 'pilotes');

% ----- Tracé des courbes TEB -----
figure('Position', [100 100 1000 700]);

semilogy(P.EbN0_dB, TEB_parfait, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
grid on; hold on;
semilogy(P.EbN0_dB, TEB_pilotes, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);

xlabel("E_b/N_0 (dB)", 'FontSize', 12);
ylabel("TEB", 'FontSize', 12);
legend("Canal parfait (CSI)", "Canal estimé (4 pilotes)", ...
    "Location", "southwest", 'FontSize', 11);
title("PARTIE 8 : Impact de l'estimation du canal par pilotes", 'FontSize', 13);
ylim([1e-5 1]);

% Annotation
text(15, 0.3, sprintf('Pilotes: indices %s', mat2str(P.idx_pilotes)), ...
    'FontSize', 10, 'BackgroundColor', 'white');
text(15, 0.15, 'Interpolation linéaire', 'FontSize', 10, 'BackgroundColor', 'white');

fprintf("\n===== Simulations terminées =====\n");

% ----- Visualisation de l'estimation du canal -----
figure('Position', [150 150 1200 500]);

% Exemple d'estimation pour un SNR élevé
% On refait une estimation pour visualiser
EbN0_test = 20;  % dB
EbN0_lin = 10^(EbN0_test/10);

% Génération d'un symbole OFDM avec pilotes (mapping IEEE 802.11a)
nb_test = 100;
Xf_test = zeros(P.Nfft, nb_test);
Nu_data = length(P.idx_data_80211a);  % 48 data
symboles_test = (2*randi([0 1], Nu_data*nb_test, 1)-1 + 1j*(2*randi([0 1], Nu_data*nb_test, 1)-1))/sqrt(2);
Xf_test(P.idx_data_80211a, :) = reshape(symboles_test, Nu_data, nb_test);
for kk = 1:length(P.idx_pilotes)
    Xf_test(P.idx_pilotes(kk), :) = P.symboles_pilotes(kk);
end

% Transmission
xt_test = ifft(Xf_test, P.Nfft);
xt_cp_test = [xt_test(end-P.Ncp+1:end, :); xt_test];
x_test = xt_cp_test(:);
puissance_test = mean(abs(x_test).^2);
y_canal_test = conv(x_test, h);
y_canal_test = y_canal_test(1:length(x_test));
Pne_test = (P.Nfft/P.Nu) * (puissance_test/(2*P.k)) * (1/EbN0_lin);
bruit_test = sqrt(Pne_test/2) * (randn(size(x_test)) + 1j*randn(size(x_test)));
y_test = y_canal_test + bruit_test;

% Réception
Ybloc_test = reshape(y_test, P.Nfft+P.Ncp, nb_test);
yt_test = Ybloc_test(P.Ncp+1:end, :);
Yf_test = fft(yt_test, P.Nfft, 1);

% Estimation
Yf_pilotes_test = Yf_test(P.idx_pilotes, :);
Hest_test = estimer_canal_pilotes(Yf_pilotes_test, P);

% Tracé
subplot(1,2,1);
plot(1:P.Nfft, abs(Heff), 'b-', 'LineWidth', 2); hold on;
plot(1:P.Nfft, abs(Hest_test), 'r--', 'LineWidth', 2);
plot(P.idx_pilotes, abs(Heff(P.idx_pilotes)), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
grid on;
xlabel("Indice sous-porteuse", 'FontSize', 11);
ylabel("|H[k]|", 'FontSize', 11);
legend("H parfait", "H estimé", "Pilotes", "Location", "best");
title(sprintf("|H[k]| - Eb/N0 = %d dB", EbN0_test), 'FontSize', 12);

subplot(1,2,2);
plot(1:P.Nfft, angle(Heff), 'b-', 'LineWidth', 2); hold on;
plot(1:P.Nfft, angle(Hest_test), 'r--', 'LineWidth', 2);
plot(P.idx_pilotes, angle(Heff(P.idx_pilotes)), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
grid on;
xlabel("Indice sous-porteuse", 'FontSize', 11);
ylabel("Phase(H[k]) (rad)", 'FontSize', 11);
legend("H parfait", "H estimé", "Pilotes", "Location", "best");
title(sprintf("Phase(H[k]) - Eb/N0 = %d dB", EbN0_test), 'FontSize', 12);

sgtitle("Comparaison canal réel vs estimé par pilotes", 'FontSize', 14);

% ----- Calcul de l'erreur d'estimation -----
erreur_estimation = norm(Hest_test - Heff) / norm(Heff);
fprintf("Erreur relative d'estimation du canal : %.2f%%\n", erreur_estimation*100);
