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

% ---------- Codage convolutif et poinçonnage ----------
% Polynômes générateurs (K=3, R=1/2 de base)
P.poly_oct = [5 7];      % polynômes en octal

% Matrices de poinçonnage pour différents rendements
% (1 = bit conservé, 0 = bit supprimé/poinçonné)
% Format : vecteur colonne lu cycliquement sur les 2 sorties du codeur
P.punc_1_2 = [1; 1];           % R=1/2 : pas de poinçonnage (tous conservés)
P.punc_2_3 = [1; 1; 1; 0];     % R=2/3 : pattern [1 1; 1 0] -> garde 3 bits sur 4
P.punc_3_4 = [1; 1; 1; 0; 0; 1]; % R=3/4 : pattern [1 1 0; 1 0 1] -> garde 4 bits sur 6

% ---------- Entrelacement IEEE 802.11a ----------
P.activer_entrelacement = false;  % true pour activer l'entrelacement

% ---------- Estimation canal par pilotes (PARTIE 8) ----------
% IEEE 802.11a : Mapping OFDM centré autour de DC (k=0)
% 64 sous-porteuses au total :
%   - 52 actives (48 données + 4 pilotes)
%   - 12 nulles (DC + bandes de garde)
%
% Notation centrée DC : k = -32 à +31
% En notation MATLAB (1 à 64) : indice = mod(k, 64) + 1
%   - DC (k=0) → indice 1 (nulle)
%   - k = +1 à +26 → indices 2 à 27
%   - k = -26 à -1 → indices 39 à 64
%   - Bandes de garde : k = +27 à +31 (indices 28-32) et k = -32 à -27 (indices 33-38)
%
% Pilotes IEEE 802.11a : k = -21, -7, +7, +21
%   - k = +7  → indice 8
%   - k = +21 → indice 22
%   - k = -7  → indice 58
%   - k = -21 → indice 44
P.idx_pilotes = [8, 22, 44, 58];  % pilotes en notation MATLAB

% Séquence pilote : {1, 1, 1, -1} (issue du registre à décalage, sortie 0)
P.symboles_pilotes = [1; 1; 1; -1];

% 52 sous-porteuses actives (indices MATLAB)
% Fréquences positives : k = +1 à +26 → indices 2 à 27
% Fréquences négatives : k = -26 à -1 → indices 39 à 64
P.idx_actives = [2:27, 39:64];  % 52 sous-porteuses actives

% 48 sous-porteuses de données = actives moins pilotes
P.idx_data_80211a = setdiff(P.idx_actives, P.idx_pilotes);  % 48 indices

% ---------- Codage LDPC ----------
% Utilisation d'une matrice LDPC standard (DVB-S2 ou 802.11n)
% Choix : N=648, K=324, R=1/2 (compatible QPSK sur 48 porteuses)
P.ldpc_N = 648;          % longueur du mot de code
P.ldpc_K = 324;          % longueur du bloc d'information
P.ldpc_R = P.ldpc_K / P.ldpc_N;  % rendement = 1/2
P.ldpc_nb_iter = 15;     % nombre d'itérations max pour le décodeur

% ---------- Courbe Eb/N0 ----------
P.EbN0_dB = 0:2:30;

% ---------- Echantillonnage (sujet) ----------
P.Fe = 20e6;             % 20 MHz
P.Te = 1/P.Fe;           % période d'échantillonnage

% ---------- Canal TGn-B (sujet) ----------
P.tau_ns = [0 10 20 30 50 80 110 140 170];                    % retards (ns)
P.p_dB   = [0 -5.4 -10.8 -16.2 -21.6 -27 -32.4 -37.8 -43.2];  % puissances (dB)
P.p_lin  = 10.^(P.p_dB/10);   % puissances en linéaire (pour génération canal)
P.seed_canal = 3;             % graine pour reproductibilité
P.p_lin  = 10.^(P.p_dB/10);                                   % puissances (lin)

% ---------- Reproductibilité ----------
P.seed_canal = 3;         % graine RNG pour fixer une réalisation de canal

end
