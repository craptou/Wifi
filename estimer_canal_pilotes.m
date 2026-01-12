function Hest = estimer_canal_pilotes(Yf_pilotes, P)
% Estimation du canal par les porteuses pilotes

Xp = P.symboles_pilotes;

% Estimation LS sur les pilotes + moyenne temporelle
Hp_tous = Yf_pilotes ./ Xp;
Hp = mean(Hp_tous, 2);

% Interpolation sur toutes les sous-porteuses
idx_p = P.idx_pilotes(:);

% Conversion indices MATLAB -> fréquences k (centrées DC)
idx_to_k = @(idx) mod(idx - 1 + 32, 64) - 32;
k_to_idx = @(k) mod(k, 64) + 1;

k_pilotes = idx_to_k(idx_p);
[k_pilotes_sorted, ordre] = sort(k_pilotes);
Hp_sorted = Hp(ordre);

k_all = (-32:31)';

% Interpolation linéaire
Hest_real = interp1(k_pilotes_sorted, real(Hp_sorted), k_all, 'linear', 'extrap');
Hest_imag = interp1(k_pilotes_sorted, imag(Hp_sorted), k_all, 'linear', 'extrap');

% Reconversion en indices MATLAB
Hest = zeros(P.Nfft, 1);
for kk = 1:length(k_all)
    idx = k_to_idx(k_all(kk));
    Hest(idx) = Hest_real(kk) + 1j * Hest_imag(kk);
end

end
