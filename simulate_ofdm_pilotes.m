function TEB = simulate_ofdm_pilotes(P, nb_symb_ofdm, h, Heff_parfait, estimation_canal)
% Simulation OFDM avec estimation du canal par pilotes

TEB = zeros(size(P.EbN0_dB));

idx_pilotes = P.idx_pilotes;
idx_data = P.idx_data_80211a;
symboles_pilotes = P.symboles_pilotes;
Nu_data = length(idx_data);

for ii = 1:length(P.EbN0_dB)

    % Génération des bits
    Nb = Nu_data * nb_symb_ofdm * P.k;
    bits_info = randi([0 1], Nb, 1);

    % Modulation QPSK
    symboles_data = qammod(bits_info, P.M, 'InputType', 'bit', 'UnitAveragePower', true);

    % Construction trame OFDM
    Xf = zeros(P.Nfft, nb_symb_ofdm);
    Xf(idx_data, :) = reshape(symboles_data, Nu_data, nb_symb_ofdm);
    for kk = 1:length(idx_pilotes)
        Xf(idx_pilotes(kk), :) = symboles_pilotes(kk);
    end

    % IFFT + CP
    xt = ifft(Xf, P.Nfft);
    xt_cp = [xt(end-P.Ncp+1:end, :); xt];
    x = xt_cp(:);
    puissance_tx = mean(abs(x).^2);

    % Canal TGn-B
    y_canal = conv(x, h);
    y_canal = y_canal(1:length(x));

    % Bruit AWGN
    EbN0 = 10^(P.EbN0_dB(ii)/10);
    Pne = (P.Nfft/P.Nu) * (puissance_tx/(2*P.k)) * (1/EbN0);
    bruit = sqrt(Pne/2) * (randn(size(x)) + 1j*randn(size(x)));
    y = y_canal + bruit;

    % Récepteur OFDM
    Ybloc = reshape(y, P.Nfft+P.Ncp, nb_symb_ofdm);
    yt = Ybloc(P.Ncp+1:end, :);
    Yf = fft(yt, P.Nfft, 1);

    % Estimation du canal
    if strcmp(estimation_canal, 'parfait')
        Heff = Heff_parfait;
    else
        Yf_pilotes = Yf(idx_pilotes, :);
        Heff = estimer_canal_pilotes(Yf_pilotes, P);
    end

    % Égalisation ZF
    Xhat_f = Yf ./ Heff;

    % Démodulation
    z = Xhat_f(idx_data, :);
    z = z(:);
    bits_recus = qamdemod(z, P.M, 'OutputType', 'bit', 'UnitAveragePower', true);

    % Calcul TEB
    TEB(ii) = mean(bits_recus ~= bits_info);
end

end
