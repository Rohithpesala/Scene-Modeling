Fs = 10;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 184;             % Length of signal
t = (0:L-1)*T;        % Time vector
size = 184;
A = csvread('idle.csv',1,6,[1,6,size,8]);
A(:,3) = -1*A(:,3);
Y = fft(A);
P2 = abs(Y(:,2)/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1);