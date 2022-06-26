function conn_Mat = non_linear_corr(Sig)

% define number of bins
L = 15;
%loop over all channels
conn_Mat = zeros(min(size(Sig)));

for i = 1:min(size((Sig)))
    X = Sig(:,i);
    [X_binned,Idx_bins] = discretize(X,L); % split X into L bins
    X_midpoints = mean([Idx_bins(1:end-1);Idx_bins(2:end)]);
    for j = 1:min(size((Sig)))
        Y = Sig(:,j);
        Y_mean = zeros(1,L);
        for k = 1:L
            Y_mean(k) = mean(Y(X_binned == k)); % calculate average of Y bins
        end
        % interpolate
        Y_i = interp1q(X_midpoints',Y_mean',X);
        % calculate connectivity
        conn_Mat(i,j) = (nansum(Y.^2)- nansum((Y-Y_i).^2))/(nansum(Y.^2));
    end
end
end
