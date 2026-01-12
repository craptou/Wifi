function bits_desentrelaces = desentrelaceur_80211a(bits_entree, Ncbps, Nbpsc)
% Désentrelacement IEEE 802.11a

nb_symboles = length(bits_entree) / Ncbps;
s = max(Nbpsc/2, 1);

% Recalcul des permutations directes
k = (0:Ncbps-1)';
i = (Ncbps/16) * mod(k, 16) + floor(k/16);
j = s * floor(i/s) + mod(i + Ncbps - floor(16*i/Ncbps), s);

% Construction de la table inverse
idx_inverse = zeros(Ncbps, 1);
for kk = 1:Ncbps
    idx_inverse(j(kk) + 1) = kk;
end

bits_desentrelaces = zeros(size(bits_entree));
for n = 1:nb_symboles
    debut = (n-1)*Ncbps + 1;
    fin = n*Ncbps;
    bloc_entree = bits_entree(debut:fin);
    bloc_sortie = zeros(Ncbps, 1);
    for jj = 1:Ncbps
        bloc_sortie(idx_inverse(jj)) = bloc_entree(jj);
    end
    bits_desentrelaces(debut:fin) = bloc_sortie;
end

end
