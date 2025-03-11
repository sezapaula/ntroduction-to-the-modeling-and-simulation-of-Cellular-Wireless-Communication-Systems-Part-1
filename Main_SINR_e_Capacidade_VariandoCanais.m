%% Limpeza de ambiente  e variáveis
clear all;
clc;
close all;
%% Solicitação interativa do número de APs
nmrAP = input('Digite o número de APs:\n (ATENÇÃO!! Apenas quadrados perfeitos)');  % Solicita ao usuário que insira o número de APs

%% Inicialização de figuras para os graficos
figure(1); %Para a CDF da SINR
hold on;
grid on;

figure(2); % Para a CDF da capacidade
hold on;
grid on;

%% Variaveis de controle 
nmrRep = 5; % Aqui, controlamos a simul. de acordo com o número de canais, ou seja, em cada conjunto de repeticao, avaliamos diferentes quantidades de canais
nChannel = 0; %Variavel de inicializacao da quantidade de canais

for g = 1:nmrRep
    %% Definição dos parâmetros gerais da simulação
    tSQUARE = 1000;  % tamanho do grid em Metros
    bT = 10^8;  % largura de banda em Hz
    nChannel = nChannel + 1;  % numero de canais
    UEpot = 1;  % potencia de transmissão do AP (em W)
    c = 10^-4;  % constante do modelo de propagaçao
    nmrUE = 14;  % numero de usuários
    k0 = 10^-20; % constante do ruído
    a = 4;  % expoente de pathloss
    
    %% Lógica da formação das posicões dos APs e usários
    APperdim = sqrt(nmrAP); % Essa linha coloca a quantidade de APs ao longo da 'linha' do Grid
    
    % Essa linha cria o vetor APcellular
    APcellular = linspace(tSQUARE / APperdim, tSQUARE, round(APperdim)) - tSQUARE / (2 * APperdim);
    
    % Essa linha forma a matriz APcellular
    APcellular = (repmat(APcellular, round(APperdim), 1) + 1j * repmat(APcellular.', 1, round(APperdim))) * 1;
    
    nmrSetups = 100;  % número de repetições da simulação
    
    % Geração de posições aleatórias para os usuários
    UElocais = (rand(nmrSetups, nmrUE) + 1i * rand(nmrSetups, nmrUE)) * tSQUARE;
    
    %% Outros parâmetros importantes
    % Cálculo de BC (largura de banda de cada canal)
    bc = bT / nChannel;
    
    % Potência do ruído
    pN = k0 * bc;
    
    % Função local para o cálculo da potência recebida
    PRecebida = @(hor_distancias)  UEpot .* (c ./ hor_distancias.^a);  % potência recebida
    
    %% Logica do calculo de SINR e alocação de usuários em canais aleatórios
    sinr = []; %Armazenamento de valores
    
    % Cálculo das potencias recebidas e alocação de canais
    for i = 1:nmrSetups
        pot_valores = zeros(nmrAP, nmrUE);  % NmrAP X NmrUE matriz que armazena a potencia recebida
        for ch = 1:nChannel
            for j = 1:nmrUE
                distancias = abs(UElocais(i, j) - APcellular(:));  % Distancia entre o usuario e todos os APs
                pot_valores(:, j) = PRecebida(distancias);  % Potencia recebida
            end
        end
        % Calculo das maiores potencias de cada AP para cada usuario
        maiores_valores = max(pot_valores, [], 1);  % Selecao da maior potencia recebida
        
        % vetor de alocação aleatória de canais para os usuarios
        usuario_chanais = randi([1, nChannel], 1, nmrUE);
        
        % Contadores para usuarios por canal
        for f = 1:nmrAP
            for u = 1:nmrUE
                % Aqui, determinamos o canal que o usuario esta usando
                canal = usuario_chanais(u);
                
                % Calcula a interferência do canal (potência de outros usuarios no mesmo canal)
                interferencia_potencia = 0;
                for ch1 = 1:nmrUE
                    if usuario_chanais(ch1) == canal && ch1 ~= u  % Condicao importante, para que nao calcule a potencia do proprio usuario como interferencia
                        interferencia_potencia = interferencia_potencia + pot_valores(f, ch1);
                    end
                end
                
                % Calculo da SINR (Interferencia + Ruido)
                sinr_value = maiores_valores(u) / (interferencia_potencia + pN);
                sinr = [sinr, sinr_value];
            end
        end
    end
    
    conv_sinr = reshape(sinr, [nmrAP * nmrSetups, nmrUE]); %Organizamos todas as SINR calculadas
    Shannon = zeros(nmrSetups, nmrUE); % Matriz que armazena a capacidade do canal
    
    for x = 1:nmrUE
        for y = 1:nmrSetups * nmrAP
            Shannon(y, x) = bc * log2(1 + conv_sinr(y, x));  %Equacao de Shannon
        end
    end
    
    Shannon_ord = reshape(Shannon, 1, []); %Vetor que organiza as capacidades calculadas
    
    % Ordena os valores em 1 coluna de SINR e Capacidade
    Snr_1c = sinr(:);
    Shannon_1c = Shannon(:);
    
    % Ordena os valores em ordem ascendente, nessessario para o calculo das
    % cdfs
    sorted_snr = sort(Snr_1c);
    sorted_capacidade = sort(Shannon_1c * 10^(-6));
    capacidade_values=linspace(0,1,length(sorted_capacidade));
    
    % Calcula as probabilidades cumulativa da SNR
    cdf_values = linspace(0, 1, length(sorted_snr));
    valor_db = pow2db(sorted_snr); %passa os valores de sinr para a escala em dB

    % Calcula os percentis gerais para todos os usuarios da repeticao
    percentil_10 = prctile(valor_db(:) , 10); 
    percentil_50 = prctile(valor_db(:) , 50);
    percentil_90 = prctile(valor_db(:) , 90);

    % Calcula os percentis gerais para todos os usuários da repetição
    percentil_10C = prctile(Shannon_1c(:) * 10^(-6), 10);  % Convertendo para Mbps
    percentil_50C = prctile(Shannon_1c(:) * 10^(-6), 50);
    percentil_90C = prctile(Shannon_1c(:) * 10^(-6), 90);


    % Exibe os percentis para cada repeticao
    fprintf('Número de canais %d:\n', g);
    fprintf('SINR 10th percentil é: %.3f dB\n', percentil_10);
    fprintf('SINR 50th percentil é: %.3f dB\n', percentil_50);
    fprintf('SINR 90th percentil é: %.3f dB\n', percentil_90);

    fprintf('Capacidade 10th percentil é: %.3f Mbps\n', percentil_10C);
    fprintf('Capacidade 50th percentil é: %.3f Mbps\n', percentil_50C);
    fprintf('Capacidade 90th percentil é: %.3f Mbps\n', percentil_90C);
    
    % Plota a CDF para o valor atual da quantide de canais
    figure(1);
    plot(valor_db, cdf_values, 'DisplayName', ['Channels: ', num2str(nChannel)], 'LineWidth', 2);
    figure(2);
    plot(sorted_capacidade, capacidade_values, 'DisplayName', ['Channels: ', num2str(nChannel)], 'LineWidth', 2);
    
    
end

% Personalizacoes do grafico
figure(1);
xlabel('SINR (dB)', 'Interpreter', 'latex');
ylabel('CDF', 'Interpreter', 'latex');
legend('Location', 'best','Interpreter','latex');
title('SINR CDF for different channel numbers', 'Interpreter', 'latex');
set(gcf, 'Renderer', 'painters');

figure(2);
xlabel('Capacity (Mbps)', 'Interpreter', 'latex');
ylabel('CDF', 'Interpreter', 'latex');
legend('Location', 'best','Interpreter','latex');
title('Capacity CDF for different channel numbers', 'Interpreter', 'latex');
set(gcf, 'Renderer', 'painters');