function Hest = estimer_canal_pilotes(Yf_pilotes, P)
%ESTIMER_CANAL_PILOTES  Estimation du canal par les porteuses pilotes (PARTIE 8)
%
% Selon IEEE 802.11a, le canal est estimé sur les 4 sous-porteuses pilotes
% puis interpolé sur toutes les sous-porteuses actives.
%
% Entrées :
%   Yf_pilotes : symboles reçus sur les sous-porteuses pilotes (4 x nb_symb)
%                après FFT, AVANT égalisation
%   P          : structure de paramètres contenant :
%                - idx_pilotes : indices des 4 pilotes [8, 22, 44, 58]
%                - symboles_pilotes : séquence pilote connue [1; 1; 1; -1]
%                - Nfft : taille FFT (64)
%
% Sortie :
%   Hest : estimation du canal sur toutes les Nfft sous-porteuses (Nfft x 1)
%          moyenné sur tous les symboles OFDM
%
% Méthode :
%   1. Estimation LS sur les pilotes : Hp[k] = Yp[k] / Xp[k]
%   2. Interpolation linéaire vers les sous-porteuses de données
%   3. Moyennage temporel sur tous les symboles OFDM (canal stationnaire)

% Nombre de symboles OFDM
nb_symb = size(Yf_pilotes, 2);

% Symboles pilotes émis (connus du récepteur)
Xp = P.symboles_pilotes;  % vecteur colonne 4x1

% ----- Estimation LS sur les pilotes -----
% Hp = Yp / Xp (division élément par élément)
% Moyenne sur tous les symboles OFDM (canal supposé stationnaire)
Hp_tous = Yf_pilotes ./ Xp;  % 4 x nb_symb
Hp = mean(Hp_tous, 2);       % 4 x 1 (moyenne temporelle)

% ----- Interpolation sur toutes les sous-porteuses -----
% Les pilotes sont aux indices [8, 22, 44, 58] en notation MATLAB
% Cela correspond aux fréquences k = [+7, +21, -21, -7]
% Pour interpoler correctement, on travaille en fréquence k = -32 à +31

% Indices des pilotes
idx_p = P.idx_pilotes(:);  % [8; 22; 44; 58]

% Conversion indices MATLAB → fréquences k (centrées DC)
% idx 1 → k=0, idx 2 → k=1, ..., idx 32 → k=31, idx 33 → k=-32, ..., idx 64 → k=-1
idx_to_k = @(idx) mod(idx - 1 + 32, 64) - 32;
k_to_idx = @(k) mod(k, 64) + 1;

% Fréquences des pilotes
k_pilotes = idx_to_k(idx_p);  % [+7, +21, -21, -7]

% Trier les pilotes par fréquence croissante
[k_pilotes_sorted, ordre] = sort(k_pilotes);
Hp_sorted = Hp(ordre);

% Fréquences de toutes les sous-porteuses
k_all = (-32:31)';

% Interpolation linéaire en fréquence
Hest_real = interp1(k_pilotes_sorted, real(Hp_sorted), k_all, 'linear', 'extrap');
Hest_imag = interp1(k_pilotes_sorted, imag(Hp_sorted), k_all, 'linear', 'extrap');

% Reconversion en indices MATLAB
Hest = zeros(P.Nfft, 1);
for kk = 1:length(k_all)
    idx = k_to_idx(k_all(kk));
    Hest(idx) = Hest_real(kk) + 1j * Hest_imag(kk);
end

end
