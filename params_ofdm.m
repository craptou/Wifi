function P = params_ofdm()
%PARAMS_OFDM  Paramètres communs du projet.

P.N    = 64;
P.Nu   = 48;
P.Ncp  = 16;

P.M    = 4;                 % QPSK
P.k    = log2(P.M);

P.EbN0dB = 0:2:30;

%Fe = 20 MHz => Ts = 50 ns
P.Fe = 20e6;
P.Ts = 1/P.Fe;

% Version simple (AWGN): ici tu mets les 48 premières à 0..Nu-1
% (Plus tard, tu remplaceras par les indices 802.11a de la norme)
P.dataIdx = 1:P.Nu;

% TGn-B (valeurs sujet)
P.tau_ns = [0 10 20 30 50 80 110 140 170];
P.p_db = [0 -5.4 -10.8 -16.2 -21.6 -27 -32.4 -37.8 -43.2];
P.p_lin = 10.^(P.p_db/10);

% Canal fixe par défaut pour comparaisons
P.fixChannel = true;
P.seedChannel = 1;

end
