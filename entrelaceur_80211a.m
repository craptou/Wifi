function bits_entrelaces = entrelaceur_80211a(bits_entree, Ncbps, Nbpsc)
% Entrelacement IEEE 802.11a

nb_symboles = length(bits_entree) / Ncbps;
s = max(Nbpsc/2, 1);

% Indices de permutation
k = (0:Ncbps-1)';
i = (Ncbps/16) * mod(k, 16) + floor(k/16);
j = s * floor(i/s) + mod(i + Ncbps - floor(16*i/Ncbps), s);
idx_sortie = j + 1;  

bits_entrelaces = zeros(size(bits_entree));
for n = 1:nb_symboles
    debut = (n-1)*Ncbps + 1;
    fin = n*Ncbps;
    bloc_entree = bits_entree(debut:fin);
    bloc_sortie = zeros(Ncbps, 1);
    for kk = 1:Ncbps
        bloc_sortie(idx_sortie(kk)) = bloc_entree(kk);
    end
    bits_entrelaces(debut:fin) = bloc_sortie;
end

end
