%% Numerical Tours - Audio Processing 
%%Sound Processing with Short Time Fourier Transform
%%Aluno: Thiago Carvalho Miranda / Matrícula: 10213079
%% Local Fourier Analysis of Sound
% Carregando um som, usando sub-sampling.
n = 1024 * 16; 
options.n = n;
[x, fs] = load_sound('bird', n);

% Tocando o som
sound(x(:)', fs);

%% Plot do sinal
clf;
plot(1:n, x);figure;
axis('tight');
set_graphic_sizes([], 20);
title('Signal');

% Dando zoom onde há grande taxa de oscilação
p = 512;
t = 1:n;
sel = n/4 + (0:p-1);
subplot(2, 1, 1);
plot(t(sel), x(sel)); axis tight;
sel = n/2 + (0:p-1);
subplot(2, 1, 2);
plot(t(sel), x(sel)); axis tight;

%% Windowing function
t = linspace(-10, 10, 2048);
eta = 1e-5;
vmin = -2;

h = double(abs(t)<1);
hf = fftshift(abs(fft(h)));
hf = log10(eta+hf);
hf = hf/max(hf);
subplot(2, 1, 1);
title('Block window');
plot(t, h); axis([-2 2, -.1, 1.1]);
subplot(2, 1, 2);
plot(t, hf); axis([-2 2, vmin, 1.1]);
title('Fourier transform');

%% Hamming Window
h = cos(t*pi()/2) .* double(abs(t) < 1);
hf = fftshift(abs(fft(h)));
hf = log10(eta+hf);
hf = hf/max(hf);
subplot(2, 1, 1);
title('Hamming window');
plot(t, h); axis([-2 2, -.1, 1.1]);
subplot(2,1,2);
plot(t, hf); axis([-2 2, vmin, 1.1]);
title('Fourier transform');

%% Haning window
h = (cos(t*pi())+1)/2 .* double(abs(t)<1);
hf = fftshift(abs(fft(h)));
hf = log10(eta+hf); hf = hf/max(hf);
subplot(2,1,1);
title('Haning window');
plot(t, h); axis([-2 2, -.1, 1.1]);
subplot(2,1,2);
plot(t, hf); axis([-2 2, vmin, 1.1]);
title('Fourier transform');

%% Haning window normalizada
h = sqrt(2)/2 * (1+cos(t*pi())) ./ sqrt( 1+cos(t*pi()).^2 ) .* double(abs(t)<1);
hf = fftshift(abs(fft(h)));
hf = log10(eta+hf); 
hf = hf/max(hf);
subplot(2,1,1);
title('Normalized Haning window');
plot(t, h); axis([-2 2, -.1, 1.1]);
subplot(2,1,2);
plot(t, hf); axis([-2 2, vmin, 1.1]);
title('Fourier transform');

%% Short time Fourier Transform
w = 64 * 2;
q = w / 2;
t = 0:3*w-1;
t1 = t-2*w;
f = w/8
g1 = sin(pi*t/w).^2 .* double(t<w);
g2 = sin(pi*t1/w).^2 .* double(t1<w & t1>=0);
g3 = g1 .* sin(t*2*pi/w*f);
g4 = g2 .* sin(t*2*pi/w*f);
subplot(2,2,1);
plot(g1); axis('tight');
title('Position 0, frequency 0');
subplot(2,2,2);
plot(g2); axis('tight');
title('Position 2*w, frequency 0');
subplot(2,2,3);
plot(g3); axis('tight');
title('Position 0, frequency w/8');
subplot(2,2,4);
plot(g4); axis('tight');
title('Position 2*w, frequency w/8');

%% Spectrograma do sinal
S = perform_stft(x, w, q, options);
imageplot(abs(S)); axis('on');figure;
plot_spectrogram(S, x);

%Medindo a energia do sinal
e = norm(x, 'fro').^2;
%Energia dos coeficientes
eS = norm(abs(S), 'fro').^2;
disp(strcat['Energy conservation (should be 1)=' num2str(e/eS)]);

%% Redução de ruído em audio
% Criando um sinal com ruído
sigma = .2;
xn = x + rand(size(x))*sigma;
% Plotando o sinal criado
subplot(2,1,1);
plot(x); axis([1 n -1.2 1.2]);
set_graphic_sizes([], 20);
title('Original signal');
subplot(2,1,2);
plot(xn); axis([1 n -1.2 1.2]);figure;
set_graphic_sizes([], 20);
title('Noisy signal');
sound(xn, fs);

Sn = perform_stft(xn,w,q, options);
SnT = perform_thresholding(Sn, 2*sigma, 'hard');

subplot(2,1,1);
plot_spectrogram(Sn);
subplot(2,1,2);
plot_spectrogram(SnT);

%% Limite de bloqueio de audio
%Removendo o ruído usando coeficientes da STFT
Sn = perform_stft(xn,w,q, options);
SnT = perform_thresholding(Sn, sigma, 'block');
subplot(2,1,1);
plot_spectrogram(Sn);
subplot(2,1,2);
plot_spectrogram(SnT);

%% Separação de fonte
clc;
clear;
n = 1024*16;
s = 3; % numero de musicas
p = 2; 
options.subsampling = 1;
x = zeros(n,3);
[x(:,1),fs] = load_sound('bird', n, options);
[x(:,2),fs] = load_sound('female', n, options);
[x(:,3),fs] = load_sound('male', n, options);

%normalizando a energia dos sinais
x = x./repmat(std(x,1), [n 1]);

%utilizando uma matriz para misturar os sinais e computando a mistura
theta = linspace(0,pi(),s+1); theta(s+1) = [];
theta(1) = .2;
M = [cos(theta); sin(theta)];
y = x*M';

% mostrando os sinais e a mistura feita
clf;
for i=1:s
    subplot(s,1,i);
    plot(x(:,i)); axis('tight');
    set_graphic_sizes([], 20);
    title(strcat('Source #',num2str(i)));
end
figure;
% mostrando a saída micro
clf;
for i=1:p
    subplot(p,1,i);
    plot(y(:,i)); axis('tight');
    set_graphic_sizes([], 20);
    title(strcat('Micro #',num2str(i)));
end
figure;

%% Análise de som local utilizando Fourier
% Setando os parâmetros da STFT
options.n = n;
w = 128;   % tamanho da janela
q = w/4;    % overlap da janela

%Computando a STFT dos sinais
clf; X = []; Y = [];
for i=1:s
    X(:,:,i) = perform_stft(x(:,i),w,q, options);
    subplot(s,1,i);
    plot_spectrogram(X(:,:,i));
    set_graphic_sizes([], 20);
    title(strcat('Source #',num2str(i)));
end
figure;

%% Estimação de direção por clustering

%Computando a posição da nuvem de pontos
mf = size(Y,1);
mt = size(Y,2);
P = reshape(Y, [mt*mf p]);
P = [real(P);imag(P)];
npts = 6000; %numero de pontos a serem mostrados

% Mostrando os pontos originais
sel = randperm(n); sel = sel(1:npts);
plot(y(sel,1), y(sel,2), '.');
axis([-1 1 -1 1]*5);
set_graphic_sizes([], 20);
title('Time domain');
figure;
% Computando os ângulos associados aos pontos e o histograma
Theta = mod(atan2(P(:,2),P(:,1)), pi());
nbins = 100;
[h,t] = hist(Theta, nbins);
h = h/sum(h);
clf;
bar(t,h); axis('tight');

%% Gabor Tight Frame Transform
% Setando o tamanho das janelas
wlist = 32*[4 8 16 32];  
L = length(wlist);
% Dando overlap na janela
K = 2;
qlist = wlist/K;
% Mostrando a redundância média
disp( strcat(['Approximate redundancy of the dictionary=' num2str(K*L) '.']) );

% Carregando um som
n = 1024*32;
options.n = n;
[x0,fs] = load_sound('glockenspiel', n);

% Computando STFT com janelas pré especificadas
options.multichannel = 0;
S = perform_stft(x0,wlist,qlist, options);
% Mostrando o spectograma
plot_spectrogram(S, x0);
%Reconstruindo os sinais usando a transformada inversa de Gabor
x1 = perform_stft(S,wlist,qlist, options);
% Checando erros
e = norm(x0-x1)/norm(x0);
disp(strcat(['Reconstruction error (should be 0) = ' num2str(e, 3)]));

%% Gabor Tight Frame Denoising
% Adicionando ruido ao sinal
sigma = .05;
x = x0 + sigma*randn(size(x0));
%Transformada STFT
S = perform_stft(x,wlist,qlist, options);
T = sigma;
ST = perform_thresholding(S, T, 'soft');
%Reconstrução
xT = perform_stft(ST,wlist,qlist, options);
% Resultado
err = snr(x0,xT);
plot_spectrogram(ST, xT);
subplot(length(ST)+1,1,1);
title(strcat(['Denoised, SNR=' num2str(err,3), 'dB']));

%% Sparsity to Improve Audio Separation
% Carregando 3 sons
n = 1024*32;
options.n = n;
s = 3; % numero de sons
p = 2; % numero de micros
options.subsampling = 1;
x = zeros(n,3);
[x(:,1),fs] = load_sound('bird', n, options);
[x(:,2),fs] = load_sound('male', n, options);
[x(:,3),fs] = load_sound('glockenspiel', n, options);
% Normalizando a energia
x = x./repmat(std(x,1), [n 1]);
% Computando a matriz de mistura
theta = linspace(0,pi(),s+1); theta(s+1) = [];
theta(1) = .2;
M = [cos(theta); sin(theta)];
y = x*M';
% Transformando o par estéreo usando o multi canal STFT
options.multichannel = 1;
%Checando reconstrução
y1 = perform_stft(S, wlist, qlist, options);
disp(strcat(['Reconstruction error (should be 0)=' num2str(norm(y-y1,'fro')/norm(y, 'fro')) '.' ]));
% Parametro de regularização
lambdaV = .2;
% Inicialização
y1 = y;
S1 = S;
niter = 100;
err = [];
% Iterações
for i=1:niter
    % gradiente
    r = y - y1;
    Sr = perform_stft(r, wlist, qlist, options);
    S1 = cell_add(S1, Sr);
    % multi canal
    y1 = perform_stft(S1,wlist,qlist, options);
end

% Criando a nuvem de pontos para o tight frame e para os coeficientes BP
% esparsos
P1 = []; P = [];
for i=1:length(S)
    Si = reshape( S1{i}, [size(S1{i},1)*size(S1{i},2) 2] );
    P1 = cat(1, P1,  Si);
    Si = reshape( S{i}, [size(S{i},1)*size(S{i},2) 2] );
    P = cat(1, P,  Si);
end
P = [real(P);imag(P)];
P1 = [real(P1);imag(P1)];

% Mostrando as duas nuvens de pontos
p = size(P,1);
m = 10000;
sel = randperm(p); sel = sel(1:m);
clf;
subplot(1,2,1);
plot( P(sel,1),P(sel,2), '.' );
title('Tight frame coefficients');
axis([-10 10 -10 10]);
subplot(1,2,2);
plot( P1(sel,1),P1(sel,2), '.' );
title('Basis Pursuit coefficients');
axis([-10 10 -10 10]);
figure;
% Computando os ângulos de maior energia
d  = sqrt(sum(P.^2,2));
d1 = sqrt(sum(P1.^2,2));
I = find( d>.2 );
I1 = find( d1>.2 );
Theta  = mod(atan2(P(I,2),P(I,1)), pi());
Theta1 = mod(atan2(P1(I1,2),P1(I1,1)), pi());
% Computando histogramas
nbins = 150;
[h,t] = hist(Theta, nbins);
h = h/sum(h);
[h1,t1] = hist(Theta1, nbins);
h1 = h1/sum(h1);
% Exibindo histogramas
subplot(2,1,1);
bar(t,h); axis('tight');
set_graphic_sizes([], 20);
title('Tight frame coefficients');
subplot(2,1,2);
bar(t1,h1); axis('tight');
set_graphic_sizes([], 20);
title('Sparse coefficients');
