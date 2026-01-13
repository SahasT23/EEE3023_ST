%Experiment 1: Repetition Code
clear();
clc();
K=1; %message length
N=3; %codeword length
D=N; %minimum Hamming distance
Rc=K/N; %code rate
message=zeros(1,K); %initialise k-bit message to all zeros
codeword=zeros(1,N); %initialise n-bit codeword to all zeros
snr=[0 1 2 3 4 5 6 7 8 9 10]; %range of signal-to-noise ratios in dB
frames=[1000 1000 10000 100000 1000000 1000000 1000000 1000000 1000000 1000000 1000000]; %number of codewords generated for each SNR
snr_uncoded=linspace(0,max(snr),length(snr)); %range of signal-to-noise ratios for BER of uncoded BPSK which we will compare our results with
Eb_N0=10.^(snr/10); %convert snr in dB to Eb/N0
ber=zeros(1,length(snr)); %initialise vector of bit-error rates for each snr value
G=ones(1,N); %Repetition code generator matrix
for s=1:length(snr)
errors=0; %initialise error counter at the start of each snr
sigma=sqrt(1/(2*Rc*Eb_N0(s))); %standard deviation of noise
for f=1:frames(s)
message=randi([0 1],1,K); %generate 1-bit message with random 0s and 1s
codeword=mod(message*G,2); %multiply message by generator matrix to obtain the n-bit repetition codeword
x=1-2*codeword; %BPSK signal obtained by mapping a bit 0 to +1 (in-phase carrier signal) and a bit 1 to -1 (out-of-phase carrier signal)
n=sigma*randn(1,length(x)); %generate white Gaussian noise
y=x+n; %add noise to BPSK signal
d=(y<0); %BPSK demodulation. If the received signal is negative then the demodulated bit is 1, else the demodulated bit is 0
%Majority Decoding
n1=nnz(d); %count the number of 1s in the received binary word d
n0=N-n1; %the number of 0s in the recevied binary word d
if n1>n0 %if the number of 1s is greater than the number 0s in d
decoded_message=1; %decoded message bit is 1
else
decoded_message=0; %else the decoded message bit is 0
end
if decoded_message~=message %if the decoded message bit is different to the message bit then this is a bit error
errors=errors+1; %add the bit error to the error counter
end
end
ber(s)=errors/(frames(s)*K); %calculate BER of the repetition code
end

% Find crossover point where coded BER becomes better than uncoded
ber_uncoded = 0.5*erfc(sqrt(Eb_N0));

crossover_idx = find(ber < ber_uncoded, 1, 'first'); % first index where coded is better

if ~isempty(crossover_idx)
    % Interpolate for more precise crossover estimate
    if crossover_idx > 1
        % Linear interpolation between adjacent points
        snr_low = snr(crossover_idx - 1);
        snr_high = snr(crossover_idx);
        
        diff_low = ber(crossover_idx - 1) - ber_uncoded(crossover_idx - 1);
        diff_high = ber(crossover_idx) - ber_uncoded(crossover_idx);
        
        crossover_snr = snr_low - diff_low * (snr_high - snr_low) / (diff_high - diff_low);
    else
        crossover_snr = snr(crossover_idx);
    end
    
    fprintf('\n=== Crossover Analysis ===\n');
    fprintf('Crossover SNR: %.2f dB\n', crossover_snr);
    fprintf('At SNR = %d dB: Coded BER = %.2e, Uncoded BER = %.2e\n', ...
            snr(crossover_idx), ber(crossover_idx), ber_uncoded(crossover_idx));
else
    fprintf('No crossover found in SNR range\n');
end

% Also display coding gain at specific BER
target_ber = 1e-3;
snr_uncoded_interp = interp1(log10(ber_uncoded), snr, log10(target_ber), 'linear', 'extrap');
snr_coded_interp = interp1(log10(ber), snr, log10(target_ber), 'linear', 'extrap');
coding_gain = snr_uncoded_interp - snr_coded_interp;

fprintf('Coding gain at BER = %.0e: %.2f dB\n', target_ber, coding_gain);

figure(1)
ber_uncoded=0.5*erfc(sqrt(Eb_N0)); %the BER for BPSK (see lecture note on digital modulation from EEE2009)
semilogy(snr_uncoded,ber_uncoded,snr,ber); %plot graph with logarithmic y-axis and linear x-axis
legend1=sprintf('Uncoded BPSK modulation');
legend2=sprintf('(%d,%d,%d) Repetition code',N,K,D);
legend(legend1,legend2);
xlabel('Eb/N0, dB');

p = 0.5*erfc(sqrt(Rc*Eb_N0)); % coded bit error probability
ber_theory = 3*p.^2.*(1-p) + p.^3; % probability of 2 or 3 errors

% Modify plot to include theoretical curve:
semilogy(snr_uncoded, ber_uncoded, snr, ber, 'o-', snr, ber_theory, '--k');
legend('Uncoded BPSK', '(3,1,3) Simulated', '(3,1,3) Theoretical');