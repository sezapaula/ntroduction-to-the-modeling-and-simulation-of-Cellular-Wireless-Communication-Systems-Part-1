clear all;
clc;
nmrRep = 8;
%nmrAP = 0;

% Inicializa uma figura para os gráficos
figure;
hold on;
grid on;

for g = 1:nmrRep
    %% Definição dos parâmetros gerais da simulação
    tSQUARE = 1000;  % tamanho do grid em M
    bT = 10^8;  % largura de banda em Hz
    nChannel = 5;  % número de canais
    UEpot = 1;  % potência de transmissão do usuário
    c = 10^-4;  % constante do modelo de propagação
    nmrUE = 14;  % número de usuários
    nmrAP = g^2;  % número de APs
    k0 = 10^-20; % constante do ruído
    a = 4;  % expoente de pathloss
    
    %% Lógica da formação das posições dos APs e usuários
    APperdim = sqrt(nmrAP); % Essa linha coloca a quantidade de APs ao longo da 'linha' do Grid
    
    % Essa linha cria o vetor APcellular
    APcellular = linspace(tSQUARE / APperdim, tSQUARE, round(APperdim)) - tSQUARE / (2 * APperdim);
    
    % Essa linha forma a matriz APcellular
    APcellular = (repmat(APcellular, round(APperdim), 1) + 1j * repmat(APcellular.', 1, round(APperdim))) * 1;
    
    nmrSetups = 500;  % número de repetições da simulação
    
    % Geração de posições aleatórias para os usuários
    UElocations = (rand(nmrSetups, nmrUE) + 1i * rand(nmrSetups, nmrUE)) * tSQUARE;
    
    %% Outros parâmetros importantes
    % Cálculo de BC (largura de banda de cada canal)
    bc = bT / nChannel;
    
    % Potência do ruído
    pN = k0 * bc;
    
    % Função para cálculo da potência recebida
    PReceiver = @(hor_distances)  UEpot .* (c ./ hor_distances.^a);  % potência recebida
    
    %% Lógica do calculo de SINR e alocação de usuários em canais aleatórios
    sinr = []; %Armazenamento de valores
    
    % Cálculo das potências recebidas e alocação de canais
    for i = 1:nmrSetups
        pot_values = zeros(nmrAP, nmrUE);  % NmrAP X NmrUE
        for ch = 1:nChannel
            for j = 1:nmrUE
                distances = abs(UElocations(i, j) - APcellular(:));  % Distância entre o usuário e todos os APs
                pot_values(:, j) = PReceiver(distances);  % Potência recebida
            end
        end
        % Cálculo das maiores potências de cada AP para cada usuário
        maiores_valores = max(pot_values, [], 1);  % Maior potência por usuário
        
        % Alocação aleatória de canais para os usuários
        user_channels = randi([1, nChannel], 1, nmrUE);
        
        % Contadores para usuários por canal
        for f = 1:nmrAP
            for u = 1:nmrUE
                % Determina o canal que o usuário 'u' está usando
                channel = user_channels(u);
                
                % Calcula a interferência do canal (potência de outros usuários no mesmo canal)
                interference_power = 0;
                for ch1 = 1:nmrUE
                    if user_channels(ch1) == channel && ch1 ~= u  % Se o canal for o mesmo, e não for o próprio usuário
                        interference_power = interference_power + pot_values(f, ch1);
                    end
                end
                
                % Cálculo da SINR (Interferência + Ruído)
                sinr_value = maiores_valores(u) / (interference_power + pN);
                sinr = [sinr, sinr_value];
            end
        end
    end
    
    conv_sinr = reshape(sinr, [nmrAP * nmrSetups, nmrUE]);
    Shannon = zeros(nmrSetups, nmrUE);
    
    for x = 1:nmrUE
        for y = 1:nmrSetups * nmrAP
            Shannon(y, x) = bc * log2(1 + conv_sinr(y, x));  % Agora `conv_sinr` é uma matriz 2D
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
    
    % Calcula as probabilidades cumulativa da SNR
    cdf_values = linspace(0, 1, length(sorted_snr));
    valor_db = pow2db(sorted_snr);

% Calcula os percentis gerais para todos os usuários da repetição
    percentile_10 = prctile(Shannon_1c(:) * 10^(-6), 10);  % Convertendo para Mbps
    percentile_50 = prctile(Shannon_1c(:) * 10^(-6), 50);
    percentile_90 = prctile(Shannon_1c(:) * 10^(-6), 90);

    % Exibe os percentis para cada repetição
    fprintf('Repetição %d:\n', g);
    fprintf('Capacidade 10th percentil é: %.3f Mbps\n', percentile_10);
    fprintf('Capacidade 50th percentil é: %.3f Mbps\n', percentile_50);
    fprintf('Capacidade 90th percentil é: %.3f Mbps\n', percentile_90);
    
    % Plota a CDF para o valor atual de nChannel com linha mais espessa
    plot(sorted_capacidade, capacidade_values, 'DisplayName', ['Numbers of APs: ', num2str(nmrAP)], 'LineWidth', 2);
    
    
end

% Personaliza o gráfico
xlabel('Mbps', 'Interpreter', 'latex');
ylabel('CDF', 'Interpreter', 'latex');
legend('Location', 'best','Interpreter','latex');
title('CDF Capacity for different number of APs', 'Interpreter', 'latex');
set(gcf, 'Renderer', 'painters');
