function P = params_ofdm()
% Paramètres OFDM IEEE 802.11a

% OFDM
P.Nfft = 64;
P.Nu = 48;
P.Ncp = 16;
P.idx_data = 1:P.Nu;

% Modulation
P.M = 4;
P.k = log2(P.M);

% Codage convolutif
P.poly_oct = [5 7];
P.punc_1_2 = [1; 1];
P.punc_2_3 = [1; 1; 1; 0];
P.punc_3_4 = [1; 1; 1; 0; 0; 1];
P.activer_entrelacement = false;

% Pilotes IEEE 802.11a
P.idx_pilotes = [8, 22, 44, 58];
P.symboles_pilotes = [1; 1; 1; -1];
P.idx_actives = [2:27, 39:64];
P.idx_data_80211a = setdiff(P.idx_actives, P.idx_pilotes);

% LDPC
P.ldpc_N = 648;
P.ldpc_K = 324;
P.ldpc_R = P.ldpc_K / P.ldpc_N;
P.ldpc_nb_iter = 15;

% SNR
P.EbN0_dB = 0:2:30;

% Echantillonnage
P.Fe = 20e6;
P.Te = 1/P.Fe;

% Canal TGn-B
P.tau_ns = [0 10 20 30 50 80 110 140 170];
P.p_dB = [0 -5.4 -10.8 -16.2 -21.6 -27 -32.4 -37.8 -43.2];
P.p_lin = 10.^(P.p_dB/10);
P.seed_canal = 3;

end
