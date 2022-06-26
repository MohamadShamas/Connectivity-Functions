function [acc,disc] = accordance_discordance(Sig)

% normalize signals
Median_Deviation = Sig-repmat(median(Sig),size(Sig,1),1);
for i  = 1:size(Sig,2)
norm_sig(i,:) = Median_Deviation(:,i)./median(abs(Median_Deviation(:,i)));
end
% derive zU and zL from normalized signal
qU = quantile(norm_sig,0.8);
qL = quantile(-norm_sig,0.8);

q_repmatU = repmat(qU,size(norm_sig,1),1);
q_repmatL = repmat(qL,size(norm_sig,1),1);

zU = norm_sig > q_repmatU;
zL = -(norm_sig < -q_repmatL);

% calculate accordnace and discordance
acc = zeros(size(norm_sig,1));
disc = zeros(size(norm_sig,1));
for i = 1:size(norm_sig,1)
    ziU = zU(:,i);
    ziL = zL(:,i);
    sigma_i = sqrt(ziU'*ziU+ziL'*ziL);
    for j = 1:min((size(norm_sig)))
        zjU = zU(:,j);
        zjL = zL(:,j);
        sigma_j = sqrt(zjU'*zjU+zjL'*zjL);
        acc(i,j) = (ziU'*zjU + ziL'*zjL)./(sigma_i*sigma_j);
        disc(i,j) = (ziU'*zjL + ziL'*zjU)./(sigma_i*sigma_j);
    end
end

end