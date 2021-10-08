clc;
% clear all;
close all;
format long;

%% Index number = 180538F
A = 5;
B = 3;
C = 8;

Ap_ = 0.03+(0.01*A);        % maximum passband ripple = 0.08 dB
Aa_ = 45+B;                 % minimum stopband attenuation = 48 dB
omiga_p1 = (C*100)+300;     % lower passband edge = 1100 rad/s (0.55pi)
omiga_p2 = (C*100)+700;     % upper passband edge = 1500 rad/s (0.75pi)
omiga_a1 = (C*100)+150;     % lower stopband edge = 950 rad/s (0.475pi)
omiga_a2 = (C*100)+800;     % upper stopband edge = 1600 rad/s (0.8pi)
omiga_s = 2*((C*100)+1200); % sampling frequency = 4000 rad/s
T = 2*pi/omiga_s;           % sampling period = 0.00157 s

Bt = min((omiga_p1-omiga_a1),(omiga_a2-omiga_p2)); % critical transition width = 100 rad/s
omiga_c1 = omiga_p1-Bt/2;                          % lower cutoff frequency = 1050 rad/s
omiga_c2 = omiga_p2+Bt/2;                          % upper cutoff frequency = 1550 rad/s

%% Choosing delta
delta_p = (10^(Ap_/20)-1/10^(Ap_/20)+1);           
delta_a = 10^(-Aa_/20);
delta = min(delta_p,delta_a);

Ap = 20*log10((1+delta)/(1-delta));     % actual passband ripple
Aa = -20*log10(delta);                  % actual stopband attenuation

%% Choosing alpha
if (Aa<=21)
    alpha = 0;
elseif ((Aa>21)&&(Aa<=50))
    alpha = 0.5842*(Aa-21)^0.4+0.07886*(Aa-21);
else
    alpha = 0.1102*(Aa-8.7);
end

%% Choosing D
if (Aa<=21)
    D = 0.9222;
else
    D = (Aa-7.95)/14.36;
end

%% Choosing order of the filter
N = ceil((omiga_s*D/Bt)+1);
if (mod(N,2)==0)    % setting N to an odd value
    N=N+1;
end
%% Generating Kaiser Window function
n = -(N-1)/2:1:(N-1)/2;                     % defining length of the filter
beta = alpha*((1-(2*(n)/(N-1)).^2).^0.5);   % generating beta
Io_alpha = bessel(alpha);                   % generating Io_alpha
Io_beta = zeros(1,N);                       % defining empty array for Io_beta

for (x = 1:N)
    Io_beta(x) = bessel(beta(x));
end

wk_nT = Io_beta/Io_alpha;   % generating Kaiser Window function
% fvtool(wk_nT);  % visualization of Kaiser Window

%% Generating the impulse response of the ideal bandpass filter
n_lower = -(N-1)/2:-1;
n_upper = 1:(N-1)/2;
h_nT_lower = (1./(n_lower*pi)).*(sin(omiga_c2*n_lower*T)-sin(omiga_c1*n_lower*T));
h_nT_upper = (1./(n_upper*pi)).*(sin(omiga_c2*n_upper*T)-sin(omiga_c1*n_upper*T));
h_0 = (2/omiga_s)*(omiga_c2-omiga_c1);
h_nT = [h_nT_lower,h_0,h_nT_upper];     % ideal bandpass filter

hw_nT = wk_nT.*h_nT;    % applying the window and generating the filter function
% fvtool(hw_nT);  % visualization of the filter

%% Generating the frequency response of the filter
[Hw_omiga, w] = freqz(hw_nT);
omiga = (w*omiga_s)/(2*pi);
log_Hw_omiga = 20*log10(abs(Hw_omiga));

%% Plotting
% figure,stem(n+floor(N/2),hw_nT);                                % causal impulse response of the filter
% title('Impulse Response');
% xlabel('n'); ylabel('Amplitude'); grid on; axis tight;
% figure,plot(omiga,log_Hw_omiga);                        % magnitude response of the filter
% title('Magnitude Response in range (0 , Ws/2)');
% xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)'); grid on; 
% figure,plot(omiga,log_Hw_omiga);                        % magnitude response of the filter in the passband
% title('Magnitude Response in Passband');
% xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)');
% xlim([omiga_c1,omiga_c2]); ylim([-0.1,0.1]); grid on;  

%% Generating the input signal
omiga_1 = omiga_a1/2;
omiga_2 = omiga_p1 + ((omiga_p2-omiga_p1)/2);
omiga_3 = omiga_a2 + (((omiga_s/2)-omiga_a2)/2);
samples = 500;
n = 0:1:samples;
x_nT = sin(omiga_1*n.*T)+sin(omiga_2*n.*T)+sin(omiga_3*n.*T);    % input signal

%% Generating the output signal by applying the filter
N_point = length(x_nT)+length(hw_nT)-1;
x_fft = fft(x_nT,N_point);
hw_fft = fft(hw_nT,N_point);
y_fft = hw_fft.*x_fft;
y_nT = ifft(y_fft,N_point);
y_nT_delayed = y_nT(floor(N/2)+1:length(y_nT)-floor(N/2));  % delayed output signal
y_nT = y_nT(1:samples+1);       % output signal for first number of samples
y_nT_expected = sin(omiga_2*n.*T);     % expected output signal (if an ideal filter is used)

%% Plotting
% figure,subplot(3,1,1);stem(n,x_nT);     % input signal
% title('Input Signal');
% xlabel('n'); ylabel('Amplitude'); grid on;
% subplot(3,1,2);stem(n,y_nT);            % output signal
% title('Output Signal');
% xlabel('n'); ylabel('Amplitude'); grid on;
% subplot(3,1,3);stem(n,y_nT_expected);   % expected output signal
% title('Expected Output Signal');
% xlabel('n'); ylabel('Amplitude'); grid on;

%% Generating the frequency response of the input signal
N_point = 2^nextpow2(2*length(x_nT));
x_fft = fft(x_nT,N_point);
x_fft = fftshift(x_fft);
X_w = T*abs(x_fft);
w_x = omiga_s*linspace(0,1,N_point)-omiga_s/2;

%% Generating the frequency response of the output signal
N_point = 2^nextpow2(2*length(y_nT_delayed));
y_fft = fft(y_nT_delayed,N_point);
y_fft = fftshift(y_fft);
Y_w = T*abs(y_fft);
w_y = omiga_s*linspace(0,1,N_point)-omiga_s/2;

%% Generating the frequency response of the expected output signal
N_point = 2^nextpow2(2*length(y_nT_expected));
y_ex_fft = fft(y_nT_expected,N_point);
y_ex_fft = fftshift(y_ex_fft);
Y_ex_w = T*abs(y_ex_fft);
w_y_ex = omiga_s*linspace(0,1,N_point)-omiga_s/2;

%% Plotting
% figure,subplot(3,1,1);plot(w_x,X_w);                    % magnitude of the frequency response of input signal
% title('Frequency Spectrum of Input Signal');
% xlabel('Frequency (rad/s)'); ylabel('Magnitude'); grid on;
% subplot(3,1,2);plot(w_y,Y_w);                           % magnitude of the frequency response of output signal
% title('Frequency Spectrum of Output Signal');
% xlabel('Frequency (rad/s)'); ylabel('Magnitude'); grid on;
% subplot(3,1,3);plot(w_y_ex,Y_ex_w);                     % magnitude of the frequency response of expected output signal
% title('Frequency Spectrum of Expected Output Signal');
% xlabel('Frequency (rad/s)'); ylabel('Magnitude'); grid on;

%% Bessel function
function Io_x = bessel(x)
    Io_x = 1;
    k = 1;
    temp = 1;
    while (temp>10^-6)                          % ignore the components that are less than 10^-6
        temp = ((1/factorial(k))*((x/2)^k))^2;
        Io_x = Io_x + temp;
        k = k+1;
    end
end