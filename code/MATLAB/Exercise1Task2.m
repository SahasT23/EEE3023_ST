%% EEE3023 Error Correcting Codes - Task 2
%  Repetition Code Comparison: N = 3, 5, 101, 1001
%  This script runs all four repetition codes and creates a combined plot

clear();
clc();

%% Define SNR range and storage arrays
snr = [0 1 2 3 4 5 6 7 8 9 10];
Eb_N0 = 10.^(snr/10);
ber_uncoded = 0.5*erfc(sqrt(Eb_N0));

% Storage for BER results of each code
ber_N3 = zeros(1, length(snr));
ber_N5 = zeros(1, length(snr));
ber_N101 = zeros(1, length(snr));
ber_N1001 = zeros(1, length(snr));

% Storage for crossover points and coding gains
crossover_points = zeros(1, 4);
coding_gains = zeros(1, 4);

%% ========================================================================
%  (3,1,3) REPETITION CODE
%  ========================================================================
fprintf('\n=== Running (3,1,3) Repetition Code ===\n');

K = 1;
N = 3;
D = N;
Rc = K/N;
G = ones(1, N);
frames = [1000 1000 10000 100000 1000000 1000000 1000000 1000000 1000000 1000000 1000000];

for s = 1:length(snr)
    errors = 0;
    sigma = sqrt(1/(2*Rc*Eb_N0(s)));
    for f = 1:frames(s)
        message = randi([0 1], 1, K);
        codeword = mod(message*G, 2);
        x = 1 - 2*codeword;
        
        n = sigma*randn(1, length(x));
        y = x + n;
        
        d = (y < 0);
        
        % Majority decoding
        n1 = nnz(d);
        n0 = N - n1;
        
        if n1 > n0
            decoded_message = 1;
        else
            decoded_message = 0;
        end
        
        if decoded_message ~= message
            errors = errors + 1;
        end
    end
    ber_N3(s) = errors/(frames(s)*K);
end

% Crossover analysis for N=3
crossover_idx = find(ber_N3 < ber_uncoded, 1, 'first');
if ~isempty(crossover_idx) && crossover_idx > 1
    snr_low = snr(crossover_idx - 1);
    snr_high = snr(crossover_idx);
    diff_low = ber_N3(crossover_idx - 1) - ber_uncoded(crossover_idx - 1);
    diff_high = ber_N3(crossover_idx) - ber_uncoded(crossover_idx);
    crossover_points(1) = snr_low - diff_low * (snr_high - snr_low) / (diff_high - diff_low);
else
    crossover_points(1) = NaN;
end

% Coding gain at BER = 1e-3
target_ber = 1e-3;
snr_uncoded_interp = interp1(log10(ber_uncoded), snr, log10(target_ber), 'linear', 'extrap');
snr_coded_interp = interp1(log10(ber_N3), snr, log10(target_ber), 'linear', 'extrap');
coding_gains(1) = snr_uncoded_interp - snr_coded_interp;

fprintf('Crossover SNR: %.2f dB\n', crossover_points(1));
fprintf('Coding gain at BER = 1e-03: %.2f dB\n', coding_gains(1));

%% ========================================================================
%  (5,1,5) REPETITION CODE
%  ========================================================================
fprintf('\n=== Running (5,1,5) Repetition Code ===\n');

K = 1;
N = 5;
D = N;
Rc = K/N;
G = ones(1, N);
frames = [1000 1000 10000 100000 1000000 1000000 1000000 1000000 1000000 1000000 1000000];

for s = 1:length(snr)
    errors = 0;
    sigma = sqrt(1/(2*Rc*Eb_N0(s)));
    for f = 1:frames(s)
        message = randi([0 1], 1, K);
        codeword = mod(message*G, 2);
        x = 1 - 2*codeword;
        
        n = sigma*randn(1, length(x));
        y = x + n;
        
        d = (y < 0);
        
        % Majority decoding
        n1 = nnz(d);
        n0 = N - n1;
        
        if n1 > n0
            decoded_message = 1;
        else
            decoded_message = 0;
        end
        
        if decoded_message ~= message
            errors = errors + 1;
        end
    end
    ber_N5(s) = errors/(frames(s)*K);
end

% Crossover analysis for N=5
crossover_idx = find(ber_N5 < ber_uncoded, 1, 'first');
if ~isempty(crossover_idx) && crossover_idx > 1
    snr_low = snr(crossover_idx - 1);
    snr_high = snr(crossover_idx);
    diff_low = ber_N5(crossover_idx - 1) - ber_uncoded(crossover_idx - 1);
    diff_high = ber_N5(crossover_idx) - ber_uncoded(crossover_idx);
    crossover_points(2) = snr_low - diff_low * (snr_high - snr_low) / (diff_high - diff_low);
else
    crossover_points(2) = NaN;
end

% Coding gain at BER = 1e-3
snr_coded_interp = interp1(log10(ber_N5), snr, log10(target_ber), 'linear', 'extrap');
coding_gains(2) = snr_uncoded_interp - snr_coded_interp;

fprintf('Crossover SNR: %.2f dB\n', crossover_points(2));
fprintf('Coding gain at BER = 1e-03: %.2f dB\n', coding_gains(2));

%% ========================================================================
%  (101,1,101) REPETITION CODE
%  ========================================================================
fprintf('\n=== Running (101,1,101) Repetition Code ===\n');

K = 1;
N = 101;
D = N;
Rc = K/N;
G = ones(1, N);
frames = [1000 1000 10000 100000 1000000 1000000 1000000 1000000 1000000 1000000 1000000];

for s = 1:length(snr)
    errors = 0;
    sigma = sqrt(1/(2*Rc*Eb_N0(s)));
    for f = 1:frames(s)
        message = randi([0 1], 1, K);
        codeword = mod(message*G, 2);
        x = 1 - 2*codeword;
        
        n = sigma*randn(1, length(x));
        y = x + n;
        
        d = (y < 0);
        
        % Majority decoding
        n1 = nnz(d);
        n0 = N - n1;
        
        if n1 > n0
            decoded_message = 1;
        else
            decoded_message = 0;
        end
        
        if decoded_message ~= message
            errors = errors + 1;
        end
    end
    ber_N101(s) = errors/(frames(s)*K);
end

% Crossover analysis for N=101
crossover_idx = find(ber_N101 < ber_uncoded, 1, 'first');
if ~isempty(crossover_idx) && crossover_idx > 1
    snr_low = snr(crossover_idx - 1);
    snr_high = snr(crossover_idx);
    diff_low = ber_N101(crossover_idx - 1) - ber_uncoded(crossover_idx - 1);
    diff_high = ber_N101(crossover_idx) - ber_uncoded(crossover_idx);
    crossover_points(3) = snr_low - diff_low * (snr_high - snr_low) / (diff_high - diff_low);
else
    crossover_points(3) = NaN;
end

% Coding gain at BER = 1e-3
snr_coded_interp = interp1(log10(ber_N101), snr, log10(target_ber), 'linear', 'extrap');
coding_gains(3) = snr_uncoded_interp - snr_coded_interp;

fprintf('Crossover SNR: %.2f dB\n', crossover_points(3));
fprintf('Coding gain at BER = 1e-03: %.2f dB\n', coding_gains(3));

%% ========================================================================
%  (1001,1,1001) REPETITION CODE
%  ========================================================================
fprintf('\n=== Running (1001,1,1001) Repetition Code ===\n');

K = 1;
N = 1001;
D = N;
Rc = K/N;
G = ones(1, N);
frames = [1000 1000 10000 100000 1000000 1000000 1000000 1000000 1000000 1000000 1000000];

for s = 1:length(snr)
    errors = 0;
    sigma = sqrt(1/(2*Rc*Eb_N0(s)));
    for f = 1:frames(s)
        message = randi([0 1], 1, K);
        codeword = mod(message*G, 2);
        x = 1 - 2*codeword;
        
        n = sigma*randn(1, length(x));
        y = x + n;
        
        d = (y < 0);
        
        % Majority decoding
        n1 = nnz(d);
        n0 = N - n1;
        
        if n1 > n0
            decoded_message = 1;
        else
            decoded_message = 0;
        end
        
        if decoded_message ~= message
            errors = errors + 1;
        end
    end
    ber_N1001(s) = errors/(frames(s)*K);
end

% Crossover analysis for N=1001
crossover_idx = find(ber_N1001 < ber_uncoded, 1, 'first');
if ~isempty(crossover_idx) && crossover_idx > 1
    snr_low = snr(crossover_idx - 1);
    snr_high = snr(crossover_idx);
    diff_low = ber_N1001(crossover_idx - 1) - ber_uncoded(crossover_idx - 1);
    diff_high = ber_N1001(crossover_idx) - ber_uncoded(crossover_idx);
    crossover_points(4) = snr_low - diff_low * (snr_high - snr_low) / (diff_high - diff_low);
else
    crossover_points(4) = NaN;
end

% Coding gain at BER = 1e-3
snr_coded_interp = interp1(log10(ber_N1001), snr, log10(target_ber), 'linear', 'extrap');
coding_gains(4) = snr_uncoded_interp - snr_coded_interp;

fprintf('Crossover SNR: %.2f dB\n', crossover_points(4));
fprintf('Coding gain at BER = 1e-03: %.2f dB\n', coding_gains(4));

%% ========================================================================
%  SUMMARY TABLE
%  ========================================================================
fprintf('\n');
fprintf('=====================================================\n');
fprintf('           REPETITION CODE COMPARISON SUMMARY        \n');
fprintf('=====================================================\n');
%% 
fprintf('Code\t\tRc\t\tEnergy Penalty\tCrossover\tCoding Gain\n');
fprintf('-----------------------------------------------------\n');
fprintf('(3,1,3)\t\t1/3\t\t%.2f dB\t\t%.2f dB\t\t%.2f dB\n', 10*log10(3), crossover_points(1), coding_gains(1));
fprintf('(5,1,5)\t\t1/5\t\t%.2f dB\t\t%.2f dB\t\t%.2f dB\n', 10*log10(5), crossover_points(2), coding_gains(2));
fprintf('(101,1,101)\t1/101\t\t%.2f dB\t\t%.2f dB\t\t%.2f dB\n', 10*log10(101), crossover_points(3), coding_gains(3));
fprintf('(1001,1,1001)\t1/1001\t\t%.2f dB\t\t%.2f dB\t\t%.2f dB\n', 10*log10(1001), crossover_points(4), coding_gains(4));
fprintf('=====================================================\n');

%% ========================================================================
%  COMBINED COMPARISON PLOT
%  ========================================================================
figure(1)
semilogy(snr, ber_uncoded, 'k-', 'LineWidth', 2); hold on;
semilogy(snr, ber_N3, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(snr, ber_N5, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(snr, ber_N101, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(snr, ber_N1001, 'm-d', 'LineWidth', 1.5, 'MarkerSize', 6);

xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('Task 2: Repetition Code Performance Comparison');
legend('Uncoded BPSK', '(3,1,3) Repetition', '(5,1,5) Repetition', ...
       '(101,1,101) Repetition', '(1001,1,1001) Repetition', ...
       'Location', 'southwest');
grid on;
axis([0 10 1e-6 1]);

%% ========================================================================
%  INDIVIDUAL PLOTS (for report figures)
%  ========================================================================
figure(2)
semilogy(snr, ber_uncoded, 'k-', 'LineWidth', 2); hold on;
semilogy(snr, ber_N3, 'b-o', 'LineWidth', 1.5);
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('(3,1,3) Repetition Code vs Uncoded BPSK');
legend('Uncoded BPSK', '(3,1,3) Repetition Code', 'Location', 'southwest');
grid on;

figure(3)
semilogy(snr, ber_uncoded, 'k-', 'LineWidth', 2); hold on;
semilogy(snr, ber_N5, 'r-s', 'LineWidth', 1.5);
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('(5,1,5) Repetition Code vs Uncoded BPSK');
legend('Uncoded BPSK', '(5,1,5) Repetition Code', 'Location', 'southwest');
grid on;

figure(4)
semilogy(snr, ber_uncoded, 'k-', 'LineWidth', 2); hold on;
semilogy(snr, ber_N101, 'g-^', 'LineWidth', 1.5);
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('(101,1,101) Repetition Code vs Uncoded BPSK');
legend('Uncoded BPSK', '(101,1,101) Repetition Code', 'Location', 'southwest');
grid on;

figure(5)
semilogy(snr, ber_uncoded, 'k-', 'LineWidth', 2); hold on;
semilogy(snr, ber_N1001, 'm-d', 'LineWidth', 1.5);
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('(1001,1,1001) Repetition Code vs Uncoded BPSK');
legend('Uncoded BPSK', '(1001,1,1001) Repetition Code', 'Location', 'southwest');
grid on;

fprintf('\nPlots generated:\n');
fprintf('  Figure 1: Combined comparison (all codes)\n');
fprintf('  Figure 2: (3,1,3) individual plot\n');
fprintf('  Figure 3: (5,1,5) individual plot\n');
fprintf('  Figure 4: (101,1,101) individual plot\n');
fprintf('  Figure 5: (1001,1,1001) individual plot\n');