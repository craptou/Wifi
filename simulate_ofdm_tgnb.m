function TEB = simulate_ofdm_tgnb(P, nb_symb_ofdm, h, Heff, type_code, rendement_code)
% Simulation OFDM sur canal TGn-B + AWGN, égalisation ZF

TEB = zeros(size(P.EbN0_dB));

% Configuration du codage
if strcmp(type_code, 'none')
    activer_codage = false;
    type_codage = 'none';
    
elseif strcmp(type_code, 'conv')
    activer_codage = true;
    type_codage = 'convolutif';
    treillis = poly2trellis(3, P.poly_oct);
    profondeur_tb = 10;
    
    if abs(rendement_code - 1/2) < 1e-6
        mat_punc = P.punc_1_2;
    elseif abs(rendement_code - 2/3) < 1e-6
        mat_punc = P.punc_2_3;
    else
        mat_punc = P.punc_3_4;
    end
    
else
    activer_codage = true;
    type_codage = 'ldpc';
    
    % Matrice LDPC R=1/2, N=648
    blockSize = 27;
    P_ldpc = [
        0 -1 -1 -1  0  0 -1 -1  0 -1 -1  0  1  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
       22  0 -1 -1 17 -1  0  0 12 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1 -1;
        6 -1  0 -1 10 -1 -1 -1 24 -1  0 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1;
        2 -1 -1  0 20 -1 -1 -1 25  0 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1;
       23 -1 -1 -1  3 -1 -1 -1  0 -1  9 11 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1;
       24 -1 23  1 17 -1  3 -1 10 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1;
       25 -1 -1 -1  8 -1 -1 -1  7 18 -1 -1  0 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1;
       13 24 -1 -1  0 -1  8 -1  6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1;
        7 20 -1 16 22 10 -1 -1 23 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1;
       11 -1 -1 -1 19 -1 -1 -1 13 -1  3 17 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1;
       25 -1  8 -1 23 18 -1 14  9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0;
        3 -1 -1 -1 16 -1 -1  2 25  5 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0
    ];
    H_ldpc = ldpcQuasiCyclicMatrix(blockSize, P_ldpc);
    cfg_ldpc_enc = ldpcEncoderConfig(H_ldpc);
    cfg_ldpc_dec = ldpcDecoderConfig(H_ldpc);
end

for ii = 1:length(P.EbN0_dB)

    % Génération des bits et codage
    if ~activer_codage
        Nb = P.Nu * nb_symb_ofdm * P.k;
        bits_info = randi([0 1], Nb, 1);
        bits_a_moduler = bits_info;
        nb_symb_utilises = nb_symb_ofdm;
        
    elseif strcmp(type_codage, 'convolutif')
        Nb_info = P.Nu * nb_symb_ofdm * P.k;
        bits_info = randi([0 1], Nb_info, 1);
        bits_codes = convenc(bits_info, treillis, mat_punc);
        
        Ncbps = P.Nu * P.k;
        nb_symb_necessaires = ceil(length(bits_codes) / Ncbps);
        Nb_cible = nb_symb_necessaires * Ncbps;
        
        if length(bits_codes) < Nb_cible
            bits_codes = [bits_codes; zeros(Nb_cible - length(bits_codes), 1)];
        end
        
        if P.activer_entrelacement
            bits_a_moduler = entrelaceur_80211a(bits_codes, Ncbps, P.k);
        else
            bits_a_moduler = bits_codes;
        end
        nb_symb_utilises = nb_symb_necessaires;
        
    else
        % LDPC
        Nb_info_total = P.Nu * nb_symb_ofdm * P.k;
        nb_blocs_ldpc = ceil(Nb_info_total / P.ldpc_K);
        Nb_info_necessaire = nb_blocs_ldpc * P.ldpc_K;
        bits_info = double(randi([0 1], Nb_info_necessaire, 1));
        Nb_info_utile = Nb_info_total;
        bits_info_originaux = bits_info;
        
        bits_codes_ldpc = zeros(nb_blocs_ldpc * P.ldpc_N, 1);
        for bloc = 1:nb_blocs_ldpc
            debut_info = (bloc-1) * P.ldpc_K + 1;
            fin_info = bloc * P.ldpc_K;
            bloc_code = ldpcEncode(double(bits_info(debut_info:fin_info)), cfg_ldpc_enc);
            debut_code = (bloc-1) * P.ldpc_N + 1;
            fin_code = bloc * P.ldpc_N;
            bits_codes_ldpc(debut_code:fin_code) = double(bloc_code);
        end
        
        Ncbps = P.Nu * P.k;
        nb_symb_necessaires = ceil(length(bits_codes_ldpc) / Ncbps);
        Nb_cible = nb_symb_necessaires * Ncbps;
        
        if length(bits_codes_ldpc) < Nb_cible
            bits_codes_ldpc = [bits_codes_ldpc; zeros(Nb_cible - length(bits_codes_ldpc), 1)];
        end
        
        if P.activer_entrelacement
            bits_a_moduler = entrelaceur_80211a(bits_codes_ldpc, Ncbps, P.k);
        else
            bits_a_moduler = bits_codes_ldpc;
        end
        nb_symb_utilises = nb_symb_necessaires;
    end

    % Modulation QPSK
    symboles = qammod(bits_a_moduler, P.M, 'InputType', 'bit', 'UnitAveragePower', true);

    % Trame OFDM
    Xf = zeros(P.Nfft, nb_symb_utilises);
    Xf(P.idx_data_80211a,:) = reshape(symboles, P.Nu, nb_symb_utilises);

    % IFFT + CP
    xt = ifft(Xf, P.Nfft);
    xt_cp = [xt(end-P.Ncp+1:end,:); xt];
    x = xt_cp(:);
    puissance_tx = mean(abs(x).^2);

    % Canal
    y_canal = conv(x, h);
    y_canal = y_canal(1:length(x));

    % Bruit AWGN
    EbN0 = 10^(P.EbN0_dB(ii)/10);
    Pne = (P.Nfft/P.Nu) * (puissance_tx/(2*P.k)) * (1/EbN0);
    bruit = sqrt(Pne/2) * (randn(size(x)) + 1j*randn(size(x)));
    y = y_canal + bruit;

    % Réception OFDM
    Ybloc = reshape(y, P.Nfft+P.Ncp, nb_symb_utilises);
    yt = Ybloc(P.Ncp+1:end,:);
    Yf = fft(yt, P.Nfft, 1);

    % Egalisation ZF
    Xhat_f = Yf ./ Heff;
    
    % Démodulation
    z = Xhat_f(P.idx_data_80211a,:);
    z = z(:);

    if strcmp(type_codage, 'ldpc')
        Heff_data = Heff(P.idx_data);
        Pne_eff_moyen = Pne / mean(abs(Heff_data).^2);
        llrs_recus = qamdemod(z, P.M, 'OutputType', 'llr', 'UnitAveragePower', true, 'NoiseVariance', Pne_eff_moyen);
        llrs_recus = double(llrs_recus);
        llrs_recus(isnan(llrs_recus)) = 0;
        llrs_recus(isinf(llrs_recus)) = sign(llrs_recus(isinf(llrs_recus))) * 30;
    else
        bits_recus = qamdemod(z, P.M, 'OutputType', 'bit', 'UnitAveragePower', true);
    end

    % Décodage
    if ~activer_codage
        bits_estimes = bits_recus;
        
    elseif strcmp(type_codage, 'convolutif')
        if P.activer_entrelacement
            Ncbps = P.Nu * P.k;
            bits_a_decoder = desentrelaceur_80211a(bits_recus, Ncbps, P.k);
        else
            bits_a_decoder = bits_recus;
        end
        bits_estimes_brut = vitdec(bits_a_decoder, treillis, profondeur_tb, 'trunc', 'hard', mat_punc);
        bits_estimes = bits_estimes_brut(1:Nb_info);
        
    else
        % LDPC
        if P.activer_entrelacement
            Ncbps = P.Nu * P.k;
            llrs_recus = desentrelaceur_80211a(llrs_recus, Ncbps, P.k);
        end
        
        bits_decodes_ldpc = zeros(nb_blocs_ldpc * P.ldpc_K, 1);
        for bloc = 1:nb_blocs_ldpc
            debut_code = (bloc-1) * P.ldpc_N + 1;
            fin_code = bloc * P.ldpc_N;
            if fin_code <= length(llrs_recus)
                llr_bloc = llrs_recus(debut_code:fin_code);
            else
                longueur_disponible = length(llrs_recus) - debut_code + 1;
                llr_bloc = [llrs_recus(debut_code:end); zeros(P.ldpc_N - longueur_disponible, 1)];
            end
            [bloc_decode, ~, ~] = ldpcDecode(double(llr_bloc), cfg_ldpc_dec, P.ldpc_nb_iter, ...
                'OutputFormat', 'info', 'DecisionType', 'hard', 'Termination', 'early');
            debut_info = (bloc-1) * P.ldpc_K + 1;
            fin_info = bloc * P.ldpc_K;
            bits_decodes_ldpc(debut_info:fin_info) = bloc_decode;
        end
        
        bits_estimes = bits_decodes_ldpc(1:Nb_info_utile);
        bits_info = bits_info_originaux(1:Nb_info_utile);
    end

    % Calcul TEB
    TEB(ii) = mean(bits_estimes ~= bits_info);
end

end
