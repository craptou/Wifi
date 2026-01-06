%%%%%%%
% LANDRY Raphaël - OFDM AWGN - QPSK - BER vs Eb/N0
%%%%%%%
clear; close all; clc;

%% ===================== PARAMS =====================
N   = 64;        % FFT size (nb sous-porteuses)
Nu  = 48;        % nb porteuses utiles (data)
Ncp = 16;        % préfixe cyclique
Nsym = 10000;      % nb symboles OFDM simulés (augmente pour lisser la courbe)

M = 4;           % QPSK
k = log2(M);     % bits / symbole QAM

EbN0dB = 0:2:30;

% Sujet: Fe = 20 MHz => Ts = 50 ns 
Fe = 20e6;
Ts = 1/Fe;

% Indices des porteuses data (version simple : Nu premières)
% (on met toutes les autres à 0, y compris pilotes, comme demandé partie AWGN)
dataIdx = 1:Nu;

% TGn-B (valeurs sujet) :contentReference[oaicite:7]{index=7}
tau_ns = [0 10 20 30 50 80 110 140 170];
p_db   = [0 -5.4 -10.8 -16.2 -21.6 -27 -32.4 -37.8 -43.2];
p_lin  = 10.^(p_db/10);

% Retards en échantillons (arrondi) + longueur h
delay_samp = round((tau_ns*1e-9)/Ts);
Lh = max(delay_samp) + 1;


%% ===================== SIMU =====================
TEB_ZF   = zeros(size(EbN0dB));
TEB_MMSE = zeros(size(EbN0dB));   % si tu veux MMSE

%% ====== Canal TGn-B FIXE pour toute la courbe (stationnaire) ======
h = zeros(Lh,1);
for kk = 1:length(tau_ns)
    beta = sqrt(p_lin(kk)/2) * (randn + 1j*randn);  % CN(0, p_lin)
    h(delay_samp(kk)+1) = h(delay_samp(kk)+1) + beta;
end

H = fft(h, N);
epsH = 1e-12;
Heff = H; 
Heff(abs(Heff) < epsH) = epsH;


for ii = 1:length(EbN0dB)

    % ----- 1) Génération bits -----
    Nb = Nu * Nsym * k;                 % bits utiles transmis
    bits = randi([0 1], Nb, 1);
    
    % ----- Codage Convalutif ----------
    trellis = poly2trellis(3, [5 7]); % Définir le treillis pour le codage
    codedBits = convenc(bits, trellis);   % Codage convalutif

    % ----- 2) Mapping QPSK (Gray, puissance unitaire) -----
    % qammod avec InputType 'bit' produit directement les symboles complexes
    symb = qammod(codedBits, M, ...
        "InputType","bit", ...
        "UnitAveragePower",true);

    % ----- 3) Grille OFDM (N x Nsym), autres porteuses à 0 -----
    Xf = zeros(N, Nsym*2);
    Xf(dataIdx,:) = reshape(symb, Nu, Nsym*2);

    % ----- 4) IFFT + CP -----
    xt = ifft(Xf, N, 1);                        % (N x Nsym)
    xt_cp = [xt(end-Ncp+1:end,:); xt];          % (N+Ncp x Nsym)
    x = xt_cp(:);                               % sérialisation

    % Puissance P du signal OFDM émis (enveloppe complexe)
    P = mean(abs(x).^2);

    %% ----- 5) Canal TGn-B : génération h[n] stationnaire (1 réalisation) -----
    % convolution + tronquage pour garder la taille
    y_chan = conv(x, h);
    y_chan = y_chan(1:length(x));

    % ----- 6) Bruit AWGN calibré Eb/N0 (formule sujet) -----
    EbN0 = 10^(EbN0dB(ii)/10);
    Pne = (N/Nu) * (P/(2*k)) * (1/EbN0);         % variance complexe (E|n|^2)
    n = sqrt(Pne/2) * (randn(size(x)) + 1j*randn(size(x)));
    y = y_chan + n;

    % ----- 7) Récepteur : reshape, suppression CP, FFT -----
    Yr = reshape(y, N+Ncp, Nsym*2);
    yt = Yr(Ncp+1:end,:);
    Yf = fft(yt, N, 1);
    
    % ----- 8) Egalisation ZF (par sous-porteuse) -----
    Xhat_ZF = Yf ./ Heff;                          % (N x Nsym)

    % ----- 9) Egalisation MMSE -----
    % snr = Es/N0 sur une sous-porteuse utile (approx). Ici Es ≈ 1 (QPSK puissance 1)
    % Pne est variance complexe du bruit sur y(t) ; 
    Wmmse = conj(Heff) ./ (abs(Heff).^2 + Pne);
    Xhat_MMSE = Wmmse .* Yf;

    % ----- 10) Extraction data + démodulation dure -----
    symb_hat_ZF   = Xhat_ZF(dataIdx,:);   
    symb_hat_ZF   = symb_hat_ZF(:);
    symb_hat_MMSE = Xhat_MMSE(dataIdx,:); 
    symb_hat_MMSE = symb_hat_MMSE(:);
    bits_hat_ZF = qamdemod(symb_hat_ZF, M, ...
        "OutputType","bit", ...
        "UnitAveragePower",true);                % décisions dures par défaut
    bits_hat_MMSE = qamdemod(symb_hat_MMSE, M, ...
        "OutputType","bit", ...
        "UnitAveragePower",true);

    % ----- Décodage de Viterbi --------
    bits_hat_ZF = vitdec(bits_hat_ZF, trellis, 34, 'trunc', 'hard');  % Décodage de Viterbi pour ZF
    bits_hat_MMSE = vitdec(bits_hat_MMSE, trellis, 34, 'trunc', 'hard'); % Décodage de Viterbi pour MMSE

    % ----- 11) TEB -----
    TEB_ZF(ii)   = mean(bits_hat_ZF   ~= bits);
    TEB_MMSE(ii) = mean(bits_hat_MMSE ~= bits);
end

% ===================== THEORIE QPSK =====================
EbN0 = 10.^(EbN0dB/10);
TEB_th = qfunc(sqrt(2*EbN0));  % QPSK Gray en AWGN

% ===================== PLOT =====================
figure;
semilogy(EbN0dB, TEB_ZF, "o-"); grid on; hold on;
semilogy(EbN0dB, TEB_MMSE, "s-");
xlabel("E_b/N_0 (dB)");
ylabel("TEB");
legend("TGn-B + ZF (canal connu)", "TGn-B + MMSE (canal connu)", "Location","southwest");
title("OFDM sur canal WiFi TGn-B (valeurs sujet) + égalisation");

figure; 
subplot(2,1,1); plot(abs(H)); grid on; title('|H[k]|');
subplot(2,1,2); plot(angle(H)); grid on; title('angle(H[k])');
