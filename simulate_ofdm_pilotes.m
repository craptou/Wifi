function TEB = simulate_ofdm_pilotes(P, nb_symb_ofdm, h, Heff_parfait, estimation_canal)
%SIMULATE_OFDM_PILOTES  Simulation OFDM avec estimation du canal par pilotes (PARTIE 8)
%
% Cette fonction simule une chaîne OFDM non codée avec :
%   - Insertion de 4 sous-porteuses pilotes
%   - Estimation du canal soit parfaite, soit par pilotes
%
% Entrées :
%   P               : paramètres (struct)
%   nb_symb_ofdm    : nombre de symboles OFDM simulés
%   h               : réponse impulsionnelle du canal (temporel)
%   Heff_parfait    : réponse fréquentielle parfaite du canal (Nfft x 1)
%   estimation_canal: 'parfait' ou 'pilotes'
%
% Sortie :
%   TEB : taux d'erreur binaire pour chaque valeur de Eb/N0

TEB = zeros(size(P.EbN0_dB));

% Indices des pilotes et des données (mapping IEEE 802.11a centré DC)
idx_pilotes = P.idx_pilotes;           % [8, 22, 44, 58]
idx_data = P.idx_data_80211a;          % 48 sous-porteuses data (excluant pilotes)
symboles_pilotes = P.symboles_pilotes; % [1; 1; 1; -1]
Nu_data = length(idx_data);            % 48

for ii = 1:length(P.EbN0_dB)

    % ============================================================
    % 1) Génération des bits (non codé)
    % ============================================================
    Nb = Nu_data * nb_symb_ofdm * P.k;
    bits_info = randi([0 1], Nb, 1);
    bits_a_moduler = bits_info;

    % ============================================================
    % 2) Modulation QPSK - puissance moyenne unitaire
    % ============================================================
    symboles_data = qammod(bits_a_moduler, P.M, ...
        'InputType', 'bit', ...
        'UnitAveragePower', true);

    % ============================================================
    % 3) Construction trame OFDM IEEE 802.11a (Nfft x nb_symb)
    %    - 48 données sur idx_data
    %    - 4 pilotes sur idx_pilotes
    %    - Reste à 0 (DC + bandes de garde)
    % ============================================================
    Xf = zeros(P.Nfft, nb_symb_ofdm);
    
    % Insertion des données sur les 48 sous-porteuses data
    symboles_data_reshape = reshape(symboles_data, Nu_data, nb_symb_ofdm);
    Xf(idx_data, :) = symboles_data_reshape;
    
    % Insertion des pilotes (mêmes symboles sur chaque symbole OFDM)
    for kk = 1:length(idx_pilotes)
        Xf(idx_pilotes(kk), :) = symboles_pilotes(kk);
    end

    % ============================================================
    % 4) IFFT + préfixe cyclique
    % ============================================================
    xt = ifft(Xf, P.Nfft);
    xt_cp = [xt(end-P.Ncp+1:end, :); xt];
    x = xt_cp(:);
    puissance_tx = mean(abs(x).^2);

    % ============================================================
    % 5) Canal TGn-B (convolution)
    % ============================================================
    y_canal = conv(x, h);
    y_canal = y_canal(1:length(x));

    % ============================================================
    % 6) Ajout bruit AWGN
    % ============================================================
    EbN0 = 10^(P.EbN0_dB(ii)/10);
    Pne = (P.Nfft/P.Nu) * (puissance_tx/(2*P.k)) * (1/EbN0);
    bruit = sqrt(Pne/2) * (randn(size(x)) + 1j*randn(size(x)));
    y = y_canal + bruit;

    % ============================================================
    % 7) Récepteur OFDM : reshape, retrait CP, FFT
    % ============================================================
    Ybloc = reshape(y, P.Nfft+P.Ncp, nb_symb_ofdm);
    yt = Ybloc(P.Ncp+1:end, :);
    Yf = fft(yt, P.Nfft, 1);

    % ============================================================
    % 8) Estimation du canal
    % ============================================================
    if strcmp(estimation_canal, 'parfait')
        % Canal parfaitement connu
        Heff = Heff_parfait;
    elseif strcmp(estimation_canal, 'pilotes')
        % Estimation par les 4 sous-porteuses pilotes
        Yf_pilotes = Yf(idx_pilotes, :);  % 4 x nb_symb
        Hest = estimer_canal_pilotes(Yf_pilotes, P);
        Heff = Hest;
    else
        error('Mode estimation_canal non reconnu : utilisez ''parfait'' ou ''pilotes''');
    end

    % ============================================================
    % 9) Égalisation ZF
    % ============================================================
    Xhat_f = Yf ./ Heff;
    
    % Protection contre NaN/Inf
    Xhat_f(isnan(Xhat_f)) = 0;
    Xhat_f(isinf(Xhat_f)) = 0;

    % ============================================================
    % 10) Démodulation (uniquement sur les données, pas les pilotes)
    % ============================================================
    z = Xhat_f(idx_data, :);
    z = z(:);

    bits_recus = qamdemod(z, P.M, ...
        'OutputType', 'bit', ...
        'UnitAveragePower', true);

    bits_estimes = bits_recus;

    % ============================================================
    % 11) Calcul TEB
    % ============================================================
    TEB(ii) = mean(bits_estimes ~= bits_info);
end

end
