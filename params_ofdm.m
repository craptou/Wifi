function P = params_ofdm()
%PARAMS_OFDM  Paramètres communs du projet.

% ---------- OFDM ----------
P.Nfft  = 64;   % taille FFT (nombre total de sous-porteuses)
P.Nu    = 48;   % nombre de sous-porteuses utiles (porteuses "data")
P.Ncp   = 16;   % longueur du préfixe cyclique

% Indices des porteuses utiles (version simplifiée : les Nu premières)
% Remarque : dans la norme 802.11a, l'allocation est centrée autour de DC
% (on fera ça plus tard).
P.idx_data = 1:P.Nu;

% ---------- Modulation ----------
P.M = 4;                 % QPSK (4-QAM)
P.k = log2(P.M);         % bits par symbole (2 en QPSK)

% ---------- Courbe Eb/N0 ----------
P.EbN0_dB = 0:2:30;

% ---------- Echantillonnage (sujet) ----------
P.Fe = 20e6;             % 20 MHz
P.Te = 1/P.Fe;           % période d'échantillonnage

% ---------- Canal TGn-B (sujet) ----------
P.tau_ns = [0 10 20 30 50 80 110 140 170];                    % retards (ns)
P.p_dB   = [0 -5.4 -10.8 -16.2 -21.6 -27 -32.4 -37.8 -43.2];  % puissances (dB)
P.p_lin  = 10.^(P.p_dB/10);                                   % puissances (lin)

% ---------- Reproductibilité ----------
P.seed_canal = 3;         % graine RNG pour fixer une réalisation de canal

end
