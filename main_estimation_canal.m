% main_estimation_canal.m - Estimation du canal par pilotes

clear; close all; clc;

P = params_ofdm();
nb_symb_ofdm = 1000;

[h, H, Heff] = tgnb_channel(P, P.seed_canal);

fprintf("Canal parfait...\n");
TEB_parfait = simulate_ofdm_pilotes(P, nb_symb_ofdm, h, Heff, 'parfait');

fprintf("Estimation par pilotes...\n");
TEB_pilotes = simulate_ofdm_pilotes(P, nb_symb_ofdm, h, Heff, 'pilotes');

% Tracé TEB
figure;
semilogy(P.EbN0_dB, TEB_parfait, 'b-o', 'LineWidth', 2); grid on; hold on;
semilogy(P.EbN0_dB, TEB_pilotes, 'r-s', 'LineWidth', 2);
xlabel("E_b/N_0 (dB)");
ylabel("TEB");
legend("Canal parfait", "Canal estimé", "Location", "southwest");
title("Impact estimation canal par pilotes");
ylim([1e-5 1]);

% Visualisation de l'estimation
figure;

EbN0_test = 20;
EbN0_lin = 10^(EbN0_test/10);

nb_test = 100;
Xf_test = zeros(P.Nfft, nb_test);
Nu_data = length(P.idx_data_80211a);
symboles_test = (2*randi([0 1], Nu_data*nb_test, 1)-1 + 1j*(2*randi([0 1], Nu_data*nb_test, 1)-1))/sqrt(2);
Xf_test(P.idx_data_80211a, :) = reshape(symboles_test, Nu_data, nb_test);
for kk = 1:length(P.idx_pilotes)
    Xf_test(P.idx_pilotes(kk), :) = P.symboles_pilotes(kk);
end

xt_test = ifft(Xf_test, P.Nfft);
xt_cp_test = [xt_test(end-P.Ncp+1:end, :); xt_test];
x_test = xt_cp_test(:);
puissance_test = mean(abs(x_test).^2);
y_canal_test = conv(x_test, h);
y_canal_test = y_canal_test(1:length(x_test));
Pne_test = (P.Nfft/P.Nu) * (puissance_test/(2*P.k)) * (1/EbN0_lin);
bruit_test = sqrt(Pne_test/2) * (randn(size(x_test)) + 1j*randn(size(x_test)));
y_test = y_canal_test + bruit_test;

Ybloc_test = reshape(y_test, P.Nfft+P.Ncp, nb_test);
yt_test = Ybloc_test(P.Ncp+1:end, :);
Yf_test = fft(yt_test, P.Nfft, 1);

Yf_pilotes_test = Yf_test(P.idx_pilotes, :);
Hest_test = estimer_canal_pilotes(Yf_pilotes_test, P);

subplot(1,2,1);
plot(1:P.Nfft, abs(Heff), 'b-'); hold on;
plot(1:P.Nfft, abs(Hest_test), 'r--');
plot(P.idx_pilotes, abs(Heff(P.idx_pilotes)), 'go', 'MarkerSize', 8);
grid on;
xlabel("Indice k");
ylabel("|H[k]|");
legend("H parfait", "H estime", "Pilotes");
title("|H[k]|");

subplot(1,2,2);
plot(1:P.Nfft, angle(Heff), 'b-'); hold on;
plot(1:P.Nfft, angle(Hest_test), 'r--');
plot(P.idx_pilotes, angle(Heff(P.idx_pilotes)), 'go', 'MarkerSize', 8);
grid on;
xlabel("Indice k");
ylabel("Phase(H[k])");
legend("H parfait", "H estime", "Pilotes");
title("Phase(H[k])");

erreur_estimation = norm(Hest_test - Heff) / norm(Heff);
fprintf("Erreur relative estimation: %.2f%%\n", erreur_estimation*100);
