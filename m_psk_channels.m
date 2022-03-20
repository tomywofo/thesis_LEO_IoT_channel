%BER teórico y simulado para 8PSK en un canal perfecto AWGN y Rayleigh 
%Es importante recalcar que el canal Rayleigh es flat fading, 
%This means that the frequency response of the channel is “flat” within the
%coherence bandwidth. En el código se asume perfect CSI, por lo que se usa
%un ecualizador zero forcing.
%La modulación, generación del canal, demodulación y detección fue hecho
%paso a paso sin utilizar objetos de MATLAB principalmente para que
%entiendan lo que hice en el código. 

%------------------------------------------------------------------------------------------------------------
%Definición de parametros de diseño
close all
clear
mkdir('figuras');
phase_offset = pi/4;
N = 10^5; % Número de bits
M = 4; %Cantidad de símbolos para 8PSK
n_bits = log2(M); %Número de bits por símbolo
No_symb = ceil(N/n_bits); %Número total de símbolos.
SNR = -2:1:30; % Relación señal ruido dB


ipHat = zeros(1,No_symb);
ipHat_awgn = zeros(1,No_symb);
nErr_im = zeros(1, length(SNR));
nErr_re = zeros(1,length(SNR));
nErr = zeros(1, length(SNR));
nErr_awgn = zeros(1, length(SNR));
nErr_awgn_re = zeros(1, length(SNR));
nErr_awgn_im = zeros(1, length(SNR));

kappa_0 = 0.5;
sigma_0 = 0.86;
s_param = 0.5;
m = -0.25;
psi_0 = 2/pi*sigma_0^2*asin(kappa_0);
f_max = 91;
f_p = 0.2*f_max;
kappa_c = 50;
N1 = 10;
N2 = 10;
if M==2
    M_str = 'B';
elseif M == 4
    M_str = 'Q';
else
    M_str = sprintf('%d-', M);
end



%----------------------------------------------------------------------------------------------------
for i = 1:length(SNR)
    SNR_montecarlo_awgn = zeros(1,21);
    SNR_montecarlo_ray = zeros(1,21);
    for j = 1:21
        % Se generan datos ya modulados
        [ip, data] = generate_PSK(No_symb, M, phase_offset);
        s = ip; % Normalización de energía de símbolo igual a 1
        n = 1/sqrt(2)*[randn(1,No_symb) + 1i*randn(1,No_symb)]; % Ruido blanco gaussiano con varianza 0
%         theta_1_n = doppler_phases(N1);
%         theta_2_n = doppler_phases(N2);
%         f_1_n = doppler_freq_jakes(f_max, N1, kappa_0);
%         f_2_n = doppler_freq_gaussian(f_max, N2, kappa_0, kappa_c);
%         [c_1_n, c_2_n] = doppler_coefficients(N1, N2, kappa_0, sigma_0);
%         time = 1:1:No_symb;
%         time = time/(36.63*f_max);
%         mu1 = sum_of_cosines(c_1_n, f_1_n, theta_1_n, time);
%         mu2 = sum_of_sines(c_1_n, f_1_n, theta_1_n, time);
%         v2 = sum_of_cosines(c_2_n, f_2_n, theta_2_n, time);
%         p = exp(s_param*v2+m);
%         theta_p = pi/2;
%         p1 = cos(2*pi*f_p*time + theta_p);
%         p2 = sin(2*pi*f_p*time + theta_p);
%         h = (p.*p1+mu1) + 1i*(p.*p2+mu2);
        h_patz = generate_patzold(kappa_0, sigma_0, s_param, m, kappa_c, N1, N2, f_max, No_symb);
        y_patz = h_patz.*s + 10^(-(SNR(i)+10*log10(n_bits))/20)*n; % Señal en el canal Rayleigh más ruido gaussiano
        y_awgn = s + 10^(-(SNR(i)+10*log10(n_bits))/20)*n; % Señal más ruido gaussiano (canal perfecto)
    %----------------------------------------------------------------------------------------------------------------

        % Ecualización perfecta. Se conoce en forma exacta el canal (perfect
        % CSI). La ecualización corresponde a cero forcing. 
        yHat = y_patz./h_patz; 
    %----------------------------------------------------------------------------------------------------------------

        % Demodulación, mediante detección de fase
        y_re = real(yHat); % Parte real
        y_im = imag(yHat); % Parte imaginaria

        phase_vec = atan2(y_im, y_re);
        % Dejar en rango 0 a 2*pi
        phase_vec = phase_vec .* (phase_vec >= 0) + (phase_vec + 2 * pi) .* (phase_vec < 0);
        
        % Se calculan fases de referencia, para encontrar la con menor
        % distancia. Se agrega 2*pi como fase de referencia, debido a la
        % discontinuidad de fase que existe cerca de 0.
        ref_phases = 2*pi*(0:M)/M+phase_offset;
        ones_matrix = ones(M+1, length(phase_vec));
        % Se genera matriz con resta de las fases de cada dato, con
        % respecto a las fases de referencia
        phase_matrix = phase_vec.*ones_matrix - ref_phases'.*ones_matrix;
        distance = abs(phase_matrix);
        % La fase que minimiza la distancia corresponde al simbolo
        [~, index] = min(distance);
        % Si la fase que minimiza era de indice 9 = (2*pi), se cambia por
        % indice 1, porque es igual a fase 0.
        index(find(index == M+1)) = 1;
        % Se pasa de codificación gray, a binario natural
        data_rec = gray2bin(index-1,'psk', M);
        % Distancia de hamming entre datos recibidos y datos enviados, para
        % calcular numero de errores
        SNR_montecarlo_ray(j) = pdist2((data_rec+1), data, 'hamming')*length(data);
    % 
    %     nErr_re(i) = size(find([real(ip)- real(ipHat)]),2); % contador de errores parte real - Rayleigh
    %     nErr_im(i) = size(find([imag(ip)- imag(ipHat)]),2); % Contador de errores parte imaginaria - Rayleigh
    %     nErr(i)=nErr_re(i)+nErr_im(i); % Total de errores Rayleigh


         % Detector para canal perfecto AWGN
        y_re_awgn = real(y_awgn); % Parte real
        y_im_awgn = imag(y_awgn); % Parte imaginaria
        %Detector para canal AWGN, es identico al anterior
        phase_vec = atan2(y_im_awgn, y_re_awgn);
        phase_vec = phase_vec .* (phase_vec >= 0) + (phase_vec + 2 * pi) .* (phase_vec < 0);
        ref_phases = 2*pi*(0:M)/M + phase_offset;
        ones_matrix = ones(M+1, length(phase_vec));
        phase_matrix = phase_vec.*ones_matrix - ref_phases'.*ones_matrix;
        distance = abs(phase_matrix);
        [vals, index] = min(distance);
        index(find(index == M+1)) = 1;
        data_rec_awgn = gray2bin(index-1,'psk', M);
        SNR_montecarlo_awgn(j) = pdist2((data_rec_awgn+1), data, 'hamming')*length(data);
    end
    if SNR(i)==-2 || SNR(i)==10 || SNR(i) == 15 || SNR(i)==30
        save_name = [sprintf('const_%spsk_', M_str), num2str(abs(SNR(i))), 'dB'];
        
%         y_patz_no_eq_im = imag(y_patz);
%         y_patz_no_eq_re = real(y_patz);
%         figure()
%         plot (y_patz_no_eq_re, y_patz_no_eq_im, '*')
%         axis([-1.5 1.5 -1.5 1.5])
%         legend('Símbolo recibido','Símbolo enviado');
%         title_str = strcat(sprintf('Constelación de %sPSK en un canal Patzold (no eq.), SNR = ', M_str), num2str(SNR(i)), ' dB');
%         title(title_str);
%         xlabel('Fase')
%         ylabel('Cuadratura')
%         saveas(gcf, ['figuras\', save_name,'_patzold_no_eq'],'png')
        
        
        figure() 
        plot(y_re,y_im,'.');
        hold on
        plot(real(s),imag(s),'r.');
        axis([-1.5 1.5 -1.5 1.5])
        legend('Símbolo recibido','Símbolo enviado');
        title_str = strcat(sprintf('Constelación de %sPSK en un canal Patzold, SNR = ', M_str), num2str(SNR(i)), ' dB');
        title(title_str);
        xlabel('Fase')
        ylabel('Cuadratura')
        saveas(gcf, ['figuras\', save_name,'_patzold'],'png')

        figure()
        plot(y_re_awgn,y_im_awgn,'.');
        hold on
        plot(real(s),imag(s),'r.');
        axis([-1.5 1.5 -1.5 1.5])
        legend('Símbolo recibido','Símbolo enviado');
        title_str = strcat(sprintf('Constelación de %sPSK en un canal AWGN, SNR = ', M_str), string(SNR(i)), ' dB');
        title(title_str)
        xlabel('Fase')
        ylabel('Cuadratura')
        saveas(gcf, ['figuras\', save_name,'_awgn'],'png')
        
    end
    % Se promedian las realizaciones de Montecarlo
    nErr(i) = mean(SNR_montecarlo_ray);
    nErr_awgn(i) = mean(SNR_montecarlo_awgn);
end

%---------------------------------------------------------------------------------------------------------------

simBer_8PSK_rayleigh = nErr/N; %BER simulado para canal Rayleigh
simBer_8PSK_awgn = nErr_awgn/N; % Ber simulado para canal AWGN
%BER teorico canal Rayleigh y canal perfercto AWGN
%theoryBer_QPSK = erf(sqrt(EbN0/2)*sin(pi/8)); % BER teórico para AWGN
theoryBer_8PSK = berawgn(SNR, 'psk', M, 'nondiff');
theoryBerRay = berfading(SNR,'psk', M, 1); % BER teórico para Rayleigh

%--------------------------------------------------------------------------------------------------------------
% Gráficas
figure
semilogy(SNR,theoryBer_8PSK,'bo-');
hold on
semilogy(SNR,theoryBerRay,'ro-');
hold on
semilogy(SNR,simBer_8PSK_awgn,'mx-');
hold on
semilogy(SNR,simBer_8PSK_rayleigh,'gx-');
axis([-2 30 10^-5 1])
grid on
legend(sprintf('Teórica %sPSK - AWGN', M_str), sprintf('Teórica %sPSK - Patzold', M_str), sprintf('Simulación %sPSK - AWGN', M_str), sprintf('Simulación %sPSK - Patzold', M_str));
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Probabilidad de error de bit para PSK')
save_name = sprintf('ber_%spsk', M_str);
saveas(gcf, ['figuras\', save_name],'png')

figure()
plot(real(s),imag(s),'r.');
axis([-1.5 1.5 -1.5 1.5])
title_str = sprintf('Símbolos antes de ser transmitidos con %sPSK', M_str);
title(title_str)
xlabel('Fase')
ylabel('Cuadratura')
saveas(gcf, ['figuras\', save_name,'_antes'],'png')
