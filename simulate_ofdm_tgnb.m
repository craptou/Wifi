function TEB_ZF = simulate_ofdm_tgnb(P, Nsym, h, Heff, useCoding, decodeType)
%SIMULATE_OFDM_TGNB  TEB sur TGn-B + AWGN avec égalisation ZF.
% Supporte: non codé, codé convolutif R=1/2 (K=3, [5 7]) en Viterbi hard/soft.
%
%   TEB_ZF = simulate_ofdm_tgnb(P, Nsym, h, Heff, useCoding, decodeType)
%
% decodeType : 'hard' ou 'soft' (uniquement si useCoding=true)

if nargin < 6 || isempty(decodeType)
    decodeType = 'hard';
end
decodeType = lower(string(decodeType));

TEB_ZF = zeros(size(P.EbN0dB));

% Code convolutif simple (K=3, R=1/2)
if useCoding
    trellis = poly2trellis(3, [5 7]);
    K = 3;
    tblen = 5*(K-1);   % 10
    if ~(decodeType=="hard" || decodeType=="soft")
        error("simulate_ofdm_tgnb:BadDecodeType", "decodeType doit etre 'hard' ou 'soft'.");
    end
end

for ii = 1:length(P.EbN0dB)

    % ===================== Bits =====================
    if ~useCoding
        Nb = P.Nu * Nsym * P.k;
        bits = randi([0 1], Nb, 1);
        bits_tx = bits;
        Nsym_use = Nsym;
    else
        % bits d'info
        Nb_info = P.Nu * Nsym * P.k;
        bits = randi([0 1], Nb_info, 1);

        % encodage R=1/2
        codedBits = convenc(bits, trellis);   % ~2*Nb_info
        bits_tx = codedBits;

        % 2x plus de bits -> 2x plus de symboles OFDM transportés
        Nsym_use = 2*Nsym;
    end

    % ===================== Modulation QPSK =====================
    symb = qammod(bits_tx, P.M, ...
        'InputType','bit', ...
        'UnitAveragePower',true);

    % ===================== Grille OFDM =====================
    Xf = zeros(P.N, Nsym_use);
    Xf(P.dataIdx,:) = reshape(symb, P.Nu, Nsym_use);

    % ===================== IFFT + CP =====================
    xt = ifft(Xf, P.N, 1);
    xt_cp = [xt(end-P.Ncp+1:end,:); xt];
    x = xt_cp(:);

    Ptx = mean(abs(x).^2);

    % ===================== Canal TGn-B =====================
    y_chan = conv(x, h);
    y_chan = y_chan(1:length(x));

    % ===================== Bruit AWGN calibré Eb/N0 =====================
    EbN0 = 10^(P.EbN0dB(ii)/10);
    Pne = (P.N/P.Nu) * (Ptx/(2*P.k)) * (1/EbN0);  % variance complexe E|n|^2
    n = sqrt(Pne/2) * (randn(size(x)) + 1j*randn(size(x)));
    y = y_chan + n;

    % ===================== RX OFDM =====================
    Yr = reshape(y, P.N+P.Ncp, Nsym_use);
    yt = Yr(P.Ncp+1:end,:);
    Yf = fft(yt, P.N, 1);

    % ===================== Egalisation ZF =====================
    Xhat = Yf ./ Heff;

    % ===================== Extraction data =====================
    z = Xhat(P.dataIdx,:); 
    z = z(:);   % symboles QPSK estimés (complexes)

    % ===================== Demod + (option) Viterbi =====================
    if ~useCoding
        % Non codé : décision dure classique
        bits_hat = qamdemod(z, P.M, ...
            'OutputType','bit', ...
            'UnitAveragePower',true);

    else
        if decodeType == "hard"
            % --- Hard: decisions 0/1 puis Viterbi hard ---
            coded_hat = qamdemod(z, P.M, ...
                'OutputType','bit', ...
                'UnitAveragePower',true);

            bits_hat = vitdec(coded_hat, trellis, tblen, 'trunc', 'hard');

        else
            nsdec = 3;                 % nb bits de quantification soft
            softMax = 2^nsdec - 1;     % valeur max
            
            % variance par dimension I/Q (approx) : Pne = E|n|^2 = E[nI^2+nQ^2]
            sigma = sqrt(Pne/2);
            
            % Observations "BPSK" sur I et Q
            rI = real(z) / sigma;
            rQ = imag(z) / sigma;
            
            % Saturation (plage typique)
            A = 4;                     % seuil de saturation (2..6 typique)
            rI = max(min(rI, A), -A);
            rQ = max(min(rQ, A), -A);
            
            % Mapping : r=+A => 0 (0 très certain), r=-A => softMax (1 très certain)
            softI = round( (A - rI) * (softMax/(2*A)) );
            softQ = round( (A - rQ) * (softMax/(2*A)) );
            
            % Mise en forme en vecteur [I Q I Q ...] attendu par vitdec
            coded_soft = zeros(2*length(z),1);
            coded_soft(1:2:end) = softI;
            coded_soft(2:2:end) = softQ;
            
            % Décodage Viterbi soft quantifié
            bits_hat = vitdec(coded_soft, trellis, tblen, 'trunc', 'soft', nsdec);
        end
    end

    % ===================== TEB =====================
    TEB_ZF(ii) = mean(bits_hat ~= bits);
end

end
