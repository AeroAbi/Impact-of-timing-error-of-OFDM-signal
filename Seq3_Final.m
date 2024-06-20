close all
clear all
N = 64;                                 % number of sub-carriers
pos = [1:15];                           % position of useful sub-carriers, central frequency = pos/T
Nu = length(pos);                       % number of useful sub-carriers
N_symb_OFDM = 1000;                     % number of symbols per sub-carrier
N_symb = Nu * N_symb_OFDM;              % total number of symbols = Nu*Nombre symbols OFDM

%channel parameters
hc = [1 0 0.5 0.25 zeros(1, N-4)];   %case1 
%hc = [1 0.5 0.25 zeros(1, N-3)];    %case2

% Define CP lengths array
CP_lengths = [1, 2, 3, 4, 5, 10];

M = 4; % QPSK
EbN0dB = [100]; % Eb/No range

% parameters for timing error analysis
%L = 20;                         % Length of channel impulse response
L=numel(hc);
M = 10;                          % Margin for synchronization error
D = L - 1 + M;                  % Length of CP (oversized CP)
N_fft = N + D;                  % Length of OFDM symbol with CP

% Generate timing errors as multiples of the sampling period
timing_errors = [-2, 0, 6];     % random timing errors
%timing_errors = [44, 69, 90];


     % Channel impulse response
    Hc = fft(hc);
    trans = conj(Hc.');
    vc = ones(1, N_symb_OFDM);
    B = trans * vc;
    
    % Equalizer coefficients
    equalizer_ML = B;
    
    % Modulation
    bits = 2 * randi(2, 1, 2 * N_symb) - 3;      % +1/-1
    symb = bits(1:2:end) + 1i * bits(2:2:end);   % mapping bits => symbols
    entree_ifft = zeros(N, N_symb_OFDM); 
    
    for ii = 1:Nu 
        entree_ifft(pos(ii) + 1, :) = symb(ii:Nu:end);    % Correct indexing for subcarriers 5-15
    end
    
     % OFDM signal, matrix form, no CP
    sortie_ifft = ifft(entree_ifft); 
    signal_ofdmnoCP = reshape(sortie_ifft, 1, N* N_symb_OFDM);

    % Plot frequency response
            figure;
            [hhfreq, w] = freqz(hc);
            plot(w/(2*pi)*N, 20*log10(abs(hhfreq)/max(abs(hhfreq))));
            grid on;
            xlabel('f/Rs');
            title('Frequency Response of the Channel');

    % Iterate over CP lengths
         for cp_length = CP_lengths
    
            % Add cyclic prefix
            CPadded = [entree_ifft(end-cp_length+1:end, :); entree_ifft];
            
            % OFDM signal, matrix form, CP
            sortie_ifft = ifft(CPadded);
            signal_ofdm = reshape(sortie_ifft, 1, (N + cp_length) * N_symb_OFDM);
        
        
       
            % Signal PSD
            [pxx, f] = pwelch(signal_ofdm(:), 1024, 512, 1024, N);
            figure
            plot(f, 10 * log10(pxx / max(abs(pxx))))
            grid on
            xlabel ('f/Rs')
            ylabel('dB')
            title (strcat('OFDM PSD , Nu= : ',num2str(Nu),' carriers', ' CP=',num2str(cp_length)))
            legend(strcat('positions = ', num2str(pos)))
            
            % Apply channel
            signal_recu = filter(hc, 1, signal_ofdm(:));
        
            % Plot PSD at receiver input
            figure;
            [pxx, f] = pwelch(signal_recu, 1024, 512, 1024, N);
            plot(f, 10*log10(pxx/max(abs(pxx))));
            grid on;
            xlabel('f/Rs');
            ylabel('dB');
            title('OFDM PSD at Receiver Input CP=',num2str(cp_length));
                
            % Demodulation
            mat_recu = reshape(signal_recu, N+cp_length, N_symb_OFDM);    % vector => matrix
            mat_demod = fft(mat_recu);                                    % FFT (demodulation) 
            mat_demod_eq = mat_demod;
            
            % Useful subcarriers recovery
            mat_demod_utile = mat_demod_eq(pos+1,:);
            
            % Apply equalization
            matrice_apres_egalisation = zeros(length(pos), N_symb_OFDM);
            for ii = 1:length(pos)
                matrice_apres_egalisation(ii, :) = mat_demod_utile(ii, :) * equalizer_ML(pos(ii));
            end
            
            figure
            hold on
            plot(real(matrice_apres_egalisation(ii, :)), imag(matrice_apres_egalisation(ii, :)), '*')
            grid on
            title(['Constellation without CP with equalization, CP=', num2str(cp_length)])    %Equalization plot without CP
            xlabel('Real part')
            ylabel('Imaginary part')
          end  

 % Plot constellations
    ii = 11;              %subcarrier 15(at position 11) for scatter diagram
    
    figure
    hold on
    plot(real(mat_demod_utile(ii, :)), imag(mat_demod_utile(ii, :)), '*')
    grid on
    title('Constellation without CP without equalizer, no noise')
    xlabel('Real part')
    ylabel('Imaginary part')

    
for case_idx = 1:3
%for case_idx = range(N_fft)
    timing_error = timing_errors(case_idx);

    % Generate CP samples
    CP_length = D;              % CP length for oversized CP or N_fft
 

    % Generate FFT window start index based on the timing error
    if case_idx == 1
        % Case 1: FFT window starts in the range [1, L-1]-FFT window starts
        % before CP
        fft_start_idx = randi([1, L-1]);
    elseif case_idx == 2
        % Case 2: FFT window starts in the range [L, D+1] -FFT starts at
        % end of CP
        % fft_start_idx = L;
        fft_start_idx = randi([L, D+1]);
    elseif case_idx == 3
        % Case 3: FFT window starts in the range [D+2, D+N]-starts at end
        % of CP adding a portion of CP to FFT window
        %fft_start_idx = D + 2;
        fft_start_idx = randi([D+2,N_fft]);
       
    end

    % Adjust FFT window start index based on timing error
    fft_start_idx = fft_start_idx + timing_error;


    % % Apply cyclic prefix 
    % ofdm_symbol_with_cp = [entree_ifft(end-CP_length+1:end), entree_ifft];
    
    

    % Apply FFT window based on the start index
    %fft_window = ofdm_symbol_with_cp(fft_start_idx : fft_start_idx + N - 1); % selects a portion of ofdm_symbol_with_cp starting from index fft_start_idx and extending to fft_start_idx + N - 1
    fft_window = symb(fft_start_idx:fft_start_idx+N-1);

    % Plot constellation without equalization
    figure;
    plot(real(fft_window(pos)), imag(fft_window(pos)), '*');
    grid on;
    xlabel('Real Part');
    ylabel('Imaginary Part');
    title(['Constellation Without Equalization - Case ', num2str(case_idx)]);
end
