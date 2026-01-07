function TEB = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, activer_codage)
%SIMULATE_OFDM_TGNB  Simulation OFDM sur canal TGn-B + AWGN, égalisation ZF.
%
% Entrées :
%   P              : paramètres (struct)
%   nb_symb_ofdm   : nombre de symboles OFDM simulés (sans codage)
%   h, Heff        : canal (fixe) en temporel et fréquentiel
%   activer_codage : true -> code convolutif + Viterbi hard
%
% Sortie :
%   TEB : vecteur TEB(Eb/N0) sur la grille P.EbN0_dB

TEB = zeros(size(P.EbN0_dB));

% ----- Codage convolutif (simple, taux 1/2) -----
if activer_codage
    treillis = poly2trellis(3, [5 7]);  % K=3, polynômes octaux (5,7)
    profondeur_tb = 5*(3-1);            % règle simple : 5*(K-1)=10
end

for ii = 1:length(P.EbN0_dB)

    % ============================================================
    % 1) Génération des bits
    % ============================================================
    if ~activer_codage
        Nb = P.Nu * nb_symb_ofdm * P.k;
        bits_info = randi([0 1], Nb, 1);
        bits_a_moduler = bits_info;
        nb_symb_utilises = nb_symb_ofdm;
    else
        Nb_info = P.Nu * nb_symb_ofdm * P.k;
        bits_info = randi([0 1], Nb_info, 1);

        bits_codes = convenc(bits_info, treillis);   % taux 1/2 => longueur doublée
        bits_a_moduler = bits_codes;

        % On transmet deux fois plus de symboles OFDM (car deux fois plus de bits)
        nb_symb_utilises = 2 * nb_symb_ofdm;
    end

    % ============================================================
    % 2) Modulation QPSK (via 4-QAM) - puissance moyenne unitaire
    % ============================================================
    symboles = qammod(bits_a_moduler, P.M, ...
        'InputType','bit', ...
        'UnitAveragePower',true);

    % ============================================================
    % 3) Construction trame OFDM en fréquence (Nfft x nb_symb)
    % ============================================================
    Xf = zeros(P.Nfft, nb_symb_utilises);
    Xf(P.idx_data,:) = reshape(symboles, P.Nu, nb_symb_utilises);

    % ============================================================
    % 4) IFFT + ajout préfixe cyclique
    % ============================================================
    xt = ifft(Xf, P.Nfft, 1);                         % (Nfft x nb_symb)
    xt_cp = [xt(end-P.Ncp+1:end,:); xt];              % (Nfft+Ncp x nb_symb)
    x = xt_cp(:);                                     % sérialisation

    puissance_tx = mean(abs(x).^2);

    % ============================================================
    % 5) Canal TGn-B (convolution)
    % ============================================================
    y_canal = conv(x, h);
    y_canal = y_canal(1:length(x));   % tronquage pour garder la même taille

    % ============================================================
    % 6) Ajout AWGN calibré via Eb/N0 (formule du sujet)
    % ============================================================
    EbN0 = 10^(P.EbN0_dB(ii)/10);

    % Pne = E{|n|^2} (variance complexe)
    Pne = (P.Nfft/P.Nu) * (puissance_tx/(2*P.k)) * (1/EbN0);

    bruit = sqrt(Pne/2) * (randn(size(x)) + 1j*randn(size(x)));
    y = y_canal + bruit;

    % ============================================================
    % 7) Récepteur OFDM : reshape, retrait CP, FFT
    % ============================================================
    Ybloc = reshape(y, P.Nfft+P.Ncp, nb_symb_utilises);
    yt = Ybloc(P.Ncp+1:end,:);
    Yf = fft(yt, P.Nfft, 1);

    % ============================================================
    % 8) Egalisation ZF par sous-porteuse
    % ============================================================
    Xhat_f = Yf ./ Heff;

    % ============================================================
    % 9) Démodulation + (option) décodage Viterbi hard
    % ============================================================
    z = Xhat_f(P.idx_data,:);
    z = z(:);

    bits_recus = qamdemod(z, P.M, ...
        'OutputType','bit', ...
        'UnitAveragePower',true);

    if activer_codage
        bits_estimes = vitdec(bits_recus, treillis, profondeur_tb, 'trunc', 'hard');
    else
        bits_estimes = bits_recus;
    end

    % ============================================================
    % 10) Calcul TEB (comparaison aux bits d'information)
    % ============================================================
    TEB(ii) = mean(bits_estimes ~= bits_info);
end

end
