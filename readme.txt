## 5 th practical 
clc;
clear all;
close all;

k = input('Enter the number of bits per symbol (e.g. 1 for BPSK, 2 for QPSK, 3 for 8PSK): ');
M = 2^k;
N = 1e5;
SNRdB = 0:2:20;

b = 0:M-1;
gray_map = bitxor(b, floor(b/2));
[~, gray_index] = sort(gray_map);
Err = zeros(1, length(SNRdB));

for i = 1:length(SNRdB)
    bits = randi([0 1], 1, N*k);
    bits_matrix = reshape(bits, k, N).';
    bin2dec_matrix = ones(N,1) * (2.^((k-1):-1:0));
    sym_dec = sum(bits_matrix .* bin2dec_matrix, 2);
    gray_sym = gray_index(sym_dec+1) - 1;
    tx_phase = 2*pi*gray_sym/M;
    s = exp(1i * tx_phase);
    noise = (randn(1,N) + 1i*randn(1,N)) / sqrt(2);
    r = s.' + 10^(-SNRdB(i)/20) * noise.';
    rx_phase = angle(r);
    rx_phase(rx_phase < 0) = rx_phase(rx_phase < 0) + 2*pi;
    rx_sym = round(rx_phase / (2*pi/M));
    rx_sym(rx_sym == M) = 0;
    bin_sym = gray_map(rx_sym + 1);
    rx_bits = dec2bin(bin_sym, k);
    rx_bits = rx_bits.';
    rx_bits = rx_bits(:).' - '0';
    Err(i) = sum(bits ~= rx_bits);
end

simBER = Err / (N*k);
theoryBER = (1/k) * erfc(sqrt(k * 10.^(SNRdB/10)) * sin(pi/M));

figure;
semilogy(SNRdB, theoryBER, 'rs-', 'LineWidth', 2);
hold on;
semilogy(SNRdB, simBER, 'bx-', 'LineWidth', 2);
grid on;
legend('Theoretical', 'Simulated');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title(['BER vs SNR for ', num2str(M), '-PSK']);
  

##scatter plot

clc; clear; close all;
M = 16; N = 1e4; SNR = 10;
x = randi([0 M-1],N,1);
y = awgn(pskmod(x,M),SNR,'measured');
scatterplot(y);



 
 ##scatterplot

clc;
clear all;
M=16;
x=[0:M-1];
N=1;
OFF=0;
z=pskmod(x,M);
figure(1)
scatterplot(z,N,OFF,"r+");
N=1;
OFF=0;
y=qammod(x,M);
figure(2)
scatterplot(y,N,OFF,"bo");


##6th 
clc;
clear all;
close all;

% --- User Input ---
choice = input('Enter 1:BPSK,2:QPSK,3:MSK,4:16-QAM,5:MPSK,6:MFSK: ');
M = 16;
N = 4; % Note: N=log2(M), hardcoded as in the photo

Eb = input('Enter Eb: ');
No = input('Enter No: ');
SNR = Eb / No; % Calculate linear SNR

% --- Probability of Error Calculation ---
Pe = 0; % Initialize Pe

if (choice < 4) % BPSK, QPSK, MSK
    Pe = 0.5 * erfc(sqrt(SNR));
end

if (choice == 4) % 16-QAM
    Pe = 2 * erfc(sqrt(0.4 * SNR));
end

if (choice == 5) % 16-MPSK
    % Using the formula from the photo
    Pe = erfc(sqrt(N * SNR)) * sin(pi / M); 
end

if (choice == 6) % 16-MFSK
    % Using the formula from the photo
    Pe = (M - 1/2) * erfc(sqrt(0.5 * N * SNR));
end

% --- Display Result ---
disp('The calculated error probability (Pe) is:');
display(Pe);


## 7th practical
clc;
clear all;
close all;
no_of_data_bits = 64;
M = 4;
n = 256;
block_size = 16;
cp_len = floor(0.1 * block_size);
data = randsrc(1, no_of_data_bits, 0:M-1);
qpsk_modulated_data = pskmod(data, M);
S2P = reshape(qpsk_modulated_data, no_of_data_bits/M, M);
number_of_subcarriers = 4;
cp_start = block_size - cp_len;
for i = 1:number_of_subcarriers
    ifft_Subcarrier(:,i) = ifft((S2P(:,i)),16);
    for j = 1:cp_len
        cyclic_prefix(j,i) = ifft_Subcarrier(j+cp_start,i);
    end
    Append_prefix(:,i) = vertcat(cyclic_prefix(:,i), ifft_Subcarrier(:,i));
end
[rows_Append_prefix, cols_Append_prefix] = size(Append_prefix);
len_ofdm_data = rows_Append_prefix * cols_Append_prefix;
ofdm_signal = reshape(Append_prefix, 1, len_ofdm_data);
channel = randn(1,2) + sqrt(-1)*randn(1,2);
after_channel = filter(channel, 1, ofdm_signal);
awgn_noise = awgn(zeros(1,length(after_channel)),0);
recvd_signal = awgn_noise + after_channel;
recvd_signal_paralleled = reshape(recvd_signal, rows_Append_prefix, cols_Append_prefix);
recvd_signal_paralleled(1:cp_len,:) = [];
for i = 1:number_of_subcarriers
    fft_data(:,i) = fft(recvd_signal_paralleled(:,i),16);
end
recvd_serial_data = reshape(fft_data, 1, (16*4));
qpsk_demodulated_data = pskdemod(recvd_serial_data,4);
figure;
stem(data);
hold on;
stem(qpsk_demodulated_data,'rx');
grid on;
xlabel('Data Points');
ylabel('Amplitude');
title('Received Signal with Error');
   


## 8 th practical


clc;
clear all;
close all;
i=input('Enter no. of elements=');
q=input('Enter joint probabilities matrix=');
sum=0;
for n=1:i
    w=0;
    for m=1:i
        p(n)=w+q(n,m);
        w=p(n);
    end
end
disp('P(x):');disp(p);
for n=1:i
    H=sum+(p(n)*log2(1/p(n)));
    sum=H;
end
disp('H(x):');disp(H);
for n=1:i
    for m=1:i
        a(n,m)=q(n,m)/p(n);
    end
end
disp('P(Y/X):');disp(a);
d=0;
for n=1:i
    for m=1:i
        if(a(n,m)>0)
            H1=d+(a(n,m)*log2(1/a(n,m)));
            d=H1;
        end
    end
end
disp('H(Y/X):');disp(H1);
m=H-H1;
disp('MI=');disp(m);
for n=1:i
    w=0;
    for m=1:i
        s(n)=w+q(m,n);
        w=s(n);
    end
end
disp('P(Y):');disp(s);
k=0;
for n=1:i
    H2=k+(s(n)*log2(1/s(n)));
    k=H2;
end
disp('H(Y):');
disp(H2);

## 9th practical 
clc;
clear all;
close all;
code_length=0;
x=input('Enter number of symbols: ');
for m=1:x
    symbols(m)=input('Enter the symbol number: ');
    p(m)=input('Enter the probability: ');
end
Hx=0;
[dict,avglen]=huffmandict(symbols,p);
hcode=huffmanenco(symbols,dict);
dsig=huffmandeco(hcode,dict);
for m=1:x
    Hx=Hx+(p(m)*(-log2(p(m))));
end
disp(Hx);
Efficiency=(Hx/avglen)*100;
disp(Efficiency);

## 10 th practical

clc;
clear;

% Get dimensions of Generator matrix G
k = input('Enter the number of message bits: ');
n = input('Enter the length of codewords: ');

fprintf('Enter the Generator matrix G (%d x %d) row-wise:\n', k, n);
G = zeros(k,n);
for i = 1:k
    G(i,:) = input(sprintf('Row %d: ', i));
end

% Get dimensions of Parity-check matrix H
r = input('Enter the number of parity-check bits (rows of H): ');
nc = input('Enter the number of columns of H (should be same as n): ');
if nc ~= n
    error('Number of columns of H must be equal to columns of G (code length n).');
end

fprintf('Enter the Parity-check matrix H (%d x %d) row-wise:\n', r, nc);
H = zeros(r,nc);
for i = 1:r
    H(i,:) = input(sprintf('Row %d: ', i));
end

% Input message vector
msg = input(sprintf('Enter the message vector of length %d: ', k));
if length(msg) ~= k
    error('Message vector length must be %d', k);
end

% Encode message: c = m * G mod 2
codeword = mod(msg * G, 2);
disp('Encoded codeword:');
disp(codeword);

% Introduce error (optional)
err_pos = input(sprintf('Enter error position to flip (1 to %d, 0 for no error): ', n));
received = codeword;
if err_pos ~= 0
    if err_pos > n || err_pos < 1
        error('Error position must be between 1 and %d or 0 for no error.', n);
    end
    received(err_pos) = mod(received(err_pos) + 1, 2); % flip bit at err_pos
end
disp('Received vector:');
disp(received);

% Compute syndrome: s = r * H' mod 2
s = mod(received * H', 2);
disp('Syndrome:');
disp(s);

% Syndrome decoding (simple error correction)
if all(s == 0)
    disp('No error detected.');
    corrected = received;
else
    % Attempt to find error position by matching syndrome in columns of H
    corrected = received;
    syndrome_found = false;
    for pos = 1:n
        if isequal(H(:,pos)', s)
            fprintf('Error detected at position %d.\n', pos);
            corrected(pos) = mod(corrected(pos) + 1, 2); % flip bit
            syndrome_found = true;
            break;
        end
    end
    if ~syndrome_found
        disp('Error pattern not found in syndrome table. Unable to correct.');
    end
end

disp('Corrected codeword:');
disp(corrected);

% Extract original message from corrected codeword (assuming systematic code)
decoded_msg = corrected(1:k);
disp('Decoded message:');
disp(decoded_msg);


10th

#########################################################################
