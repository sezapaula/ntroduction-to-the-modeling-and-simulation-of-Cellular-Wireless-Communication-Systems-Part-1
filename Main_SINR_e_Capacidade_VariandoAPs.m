%%Limpeza de ambiente e variaveis
clear all;
clc;
close all;

%% Solicitacao interativa do numero de canais
nChannel=input('Digite o número de canais:');

%% Inicializacao de figuras para os graficos
figure(1);
hold on;
grid on;

figure(2);
hold on;
grid on;

%% Variavel de controle para as repeticoes
nmrRep = 8; %Aqui, os conjuntos de repeticoes variam com a quantidade de APs

for g = 1:nmrRep
    %% Definicao dos parametros gerais da simulacao
    tSQUARE = 1000;  % tamanho do grid em M
    bT = 10^8;  % largura de banda em Hz
    UEpot = 1;  % potencia de transmissao do AP
    c = 10^-4;  % constante do modelo de propagacao
    nmrUE = 13;  % numero de usuarios
    nmrAP = g^2;  % numero de APs (apenas quadrados perfeitos)
    k0 = 10^-20; % constante do ruido
    a = 4;  % expoente de pathloss
    
    %% Logica da formacao das posicoes dos APs e usuarios
    APperdim = sqrt(nmrAP); % Essa linha coloca a quantidade de APs ao longo da 'linha' do Grid
    
    % Essa linha cria o vetor APcellular
    APcellular = linspace(tSQUARE / APperdim, tSQUARE, round(APperdim)) - tSQUARE / (2 * APperdim);
    
    % Essa linha forma a matriz APcellular
    APcellular = (repmat(APcellular, round(APperdim), 1) + 1j * repmat(APcellular.', 1, round(APperdim))) * 1;
    
    nmrSetups = 500;  % numero de repetições da simulacao
    
    % Geracao de posicoes aleatorias para os usuarios
    UElocais = (rand(nmrSetups, nmrUE) + 1i * rand(nmrSetups, nmrUE)) * tSQUARE;
    
    %% Outros parametros importantes
    % Cálculo de BC (largura de banda de cada canal)
    bc = bT / nChannel;
    
    % Potencia do ruido
    pN = k0 * bc;
    
    % Funcao para calculo da potencia recebida
    PRecebida = @(hor_distances)  UEpot .* (c ./ hor_distances.^a);  % potencia recebida
    
    %% Logica do calculo de SINR e alocação de usuarios em canais aleatorios
    sinr = []; %Armazenamento de valores
    
    % Calculo das potencias recebidas e alocacao de canais
    for i = 1:nmrSetups
        pot_valores = zeros(nmrAP, nmrUE);  % NmrAP X NmrUE
        for ch = 1:nChannel
            for j = 1:nmrUE
                distancias = abs(UElocais(i, j) - APcellular(:));  % Distancia entre o usuario e todos os APs
                pot_valores(:, j) = PRecebida(distancias);  % Potencia recebida
            end
        end
        % Calculo das maiores potencias de cada AP para cada usuario
        maiores_valores = max(pot_valores, [], 1);  % Maior potencia por usuario
        
        % Alocacao aleatoria de canais para os usuarios
        usuario_canais = randi([1, nChannel], 1, nmrUE);
        
        % Contadores para usuarios por canal
        for f = 1:nmrAP
            for u = 1:nmrUE
                % Determina o canal que o usuario esta usando
                channel = usuario_canais(u);
                
                % Calcula a interferencia do canal (potencia de outros usuarios no mesmo canal)
                interferencia_potencia = 0;
                for ch1 = 1:nmrUE
                    if usuario_canais(ch1) == channel && ch1 ~= u  % Se o canal for o mesmo, e nao for o proprio usuario
                        interferencia_potencia = interferencia_potencia + pot_valores(f, ch1);
                    end
                end
                
                % Calculo da SINR (Interferencia + Ruido)
                sinr_value = maiores_valores(u) / (interferencia_potencia + pN);
                sinr = [sinr, sinr_value];
            end
        end
    end
    
    conv_sinr = reshape(sinr, [nmrAP * nmrSetups, nmrUE]);
    Shannon = zeros(nmrSetups, nmrUE);
    
    for x = 1:nmrUE
        for y = 1:nmrSetups * nmrAP
            Shannon(y, x) = bc * log2(1 + conv_sinr(y, x));  % Equação de Shannon
        end
    end
    
    Shannon_ord = reshape(Shannon, 1, []);
    
    % Ordena os valores em 1 coluna
    Snr_1c = sinr(:);
    Shannon_1c = Shannon(:);
    
    % Ordena os valores em ordem ascendente
    sorted_snr = sort(Snr_1c);
    sorted_capacidade = sort(Shannon_1c * 10^(-6));
    capacidade_values=linspace(0,1,length(sorted_capacidade));
    
    % Calcula as probabilidades cumulativa da SINR
    cdf_values = linspace(0, 1, length(sorted_snr));
    valor_db = pow2db(sorted_snr);


     % Calcula os percentis gerais para todos os usuários da repetição
    percentil_10 = prctile(valor_db(:) , 10);  % Convertendo para Mbps
    percentil_50 = prctile(valor_db(:) , 50);
    percentil_90 = prctile(valor_db(:) , 90);

    percentil_10C = prctile(Shannon_1c(:) * 10^(-6), 10);  % Convertendo para Mbps
    percentil_50C = prctile(Shannon_1c(:) * 10^(-6), 50);
    percentil_90C = prctile(Shannon_1c(:) * 10^(-6), 90);


    % Exibe os percentis para cada repetição
    fprintf('Numero de APs: %d\n', g^2);
    fprintf('SINR 10th percentil é: %.3f dB\n', percentil_10);
    fprintf('SINR 50th percentil é: %.3f dB\n', percentil_50);
    fprintf('SINR 90th percentil é: %.3f dB\n', percentil_90);

    fprintf('Capacidade 10th percentil é: %.3f Mbps\n', percentil_10C);
    fprintf('Capacidade 50th percentil é: %.3f Mbps\n', percentil_50C);
    fprintf('Capacidade 90th percentil é: %.3f Mbps\n', percentil_90C);
   

     % Plota a CDF 
    figure(1);
    plot(valor_db, cdf_values, 'DisplayName', ['Numbers of APs ', num2str(nmrAP)], 'LineWidth', 2);
    figure(2);
    plot(sorted_capacidade, capacidade_values, 'DisplayName', ['Numbers of APs: ', num2str(nmrAP)], 'LineWidth', 2);
    
end

% Personalizacoes do grafico
figure(1);
xlabel('SINR (dB)', 'Interpreter', 'latex');
ylabel('CDF', 'Interpreter', 'latex');
legend('Location', 'best','Interpreter','latex');
title('SINR CDF for differen numbers of APs', 'Interpreter', 'latex');
set(gcf, 'Renderer', 'painters');

figure(2);
xlabel('Capacity (Mbps)', 'Interpreter', 'latex');
ylabel('CDF', 'Interpreter', 'latex');
legend('Location', 'best','Interpreter','latex');
title('Capacity CDF for different numbers of APs', 'Interpreter', 'latex');
set(gcf, 'Renderer', 'painters');