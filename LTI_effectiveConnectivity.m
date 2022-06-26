function [eff_conn_mat,MSE] = LTI_effectiveConnectivity(sigMat,step_size,Fs,plot_flag)

% [eff_conn_mat,MSE] = LTI_effectiveConnectivity(sigMat,step_size,Fs,plot_flag)
%
% Caculate effective connectivity for EEG and FMRI
%
%   INPUT
%       sigMat    = signal matrix of dimension NxS (N: number of cahnnels, S
%                   number of samples)
%       step_size = number of seconds of each window (float)
%       Fs        = sampling Frequency (Hz)
%       plot_flag = 1: plot randomly selected fitted data , 0: don't plot
%
%   OUTPUT
%       eff_conn_mat   = effective connectivity matrix A of size NxN such   
%                        that x(t+1) = A*x(t)
%       MSE            = mean square error between estimated x(t) and
%                         actual signal of size Nb_windows x N
%   Mohamad Shamas, UCLA, 2022


% tolerance value 
tol = 1e-4;
% max iteration value
max_i = 500;
% initialize number of plots
nb_plots = 5;
% get number of channels
N = size(sigMat,1);
% calculate window size T
T = floor(step_size*Fs);

% divide the signal into equal windows
nb_windows = floor(size(sigMat,2)/T);
% dismiss last window if length is not a multiple of step_size
sigMat_trunc = sigMat(:,1:nb_windows*floor(step_size*Fs));
% initiate effective connectivity matrix
eff_conn_all = zeros(nb_windows,N,N);
% loop over all windows
for win_idx = 1:nb_windows
    
    % intiate matrix H for every window
    H = zeros((T-1)*N,N*N);
    % construct window signal to be used in solving for X in b = HX
    wind_sig = sigMat_trunc(:,T*(win_idx-1)+1:T*win_idx);
    % construct b in form = [x1(2) x2(2) ... xN(2) ... x1(T) x2(T) ... xN(T)]
    b = wind_sig(:,2:end); b = b(:);
    % construct matrix H
    temp_vect = repmat(wind_sig(:,1:T-1),N,1);
    for i = 1:T-1
        H((i-1)*N+1:i*N,1:N) = reshape(temp_vect(:,i),N,[])';
    end
    circ_shift_amount = repmat((0:N-1).*N,1,T-1);
    for j = 1:N*(T-1)
        H(j,:) = circshift(H(j,:),circ_shift_amount(j));
    end
    
    % solve for X = [A11 A12 ... A1N A21 ... AN,N−1 ANN ]'
    X = lsqr(H,b,tol,max_i); 
    A = reshape(X,N,N)'; % reshape to NxN matrixl
    eff_conn_all(win_idx,:,:)= A;
    % construct estimated signal
    sig_est(:,1) = wind_sig(:,1);
    sig2plot(win_idx,:,1) = sig_est(:,1);
    for sample_idx = 1:T-1
        sig_est(:,sample_idx+1)= A*sig_est(:,sample_idx);
        sig2plot(win_idx,:,sample_idx+1) = sig_est(:,sample_idx+1);
    end
    % calculate error of estimated signal
    MSE(win_idx,:) = mean((wind_sig' - sig_est').^2);    
end
    eff_conn_mat = squeeze(mean(eff_conn_all));
    
    % plot random signals if plot_flag = true
    if plot_flag
        % choose 5 random channels
        ch_idx = randi(N,min(N,nb_plots));
        figure
        for sub_plot_idx = 1:length(ch_idx)
            subplot(length(ch_idx),1,sub_plot_idx)
            t = (0:size(sigMat_trunc,2)-1)./Fs;
            scatter(t,sigMat_trunc(ch_idx(sub_plot_idx),:),'.r');
            hold on
            sig = squeeze(sig2plot(:,ch_idx(sub_plot_idx),:))';
            sig =sig(:);
            plot(t,sig,'k');
            xlabel('time (s)')
            ylabel(['Channel ' num2str(ch_idx(sub_plot_idx))]);
        end
        
    end
end