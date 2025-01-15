%% Import e preprocessing dei dati in matlab per le prove prima del forcasting di deep
clc
clear all
data_2021 = xlsread('Dataset-Project-Deep-Learning-SMRES.xls');
data_2022 = xlsread('Dataset-Project-Deep-Learning-SMRES - 2022.xls');
data_2022 = data_2022(:,2:5);


%% calcolo della deviazione standard della potenza richiesta dagli uffici e dell'irraggiamento
devstd_2021_uffici = std(data_2021(:,1));
devstd_2022_uffici = std(data_2022(:,1));
devstd_2021_irraggiamento = std(data_2021(:,4));
devstd_2022_irraggiamento = std(data_2022(:,4));


%% aggiunta del rumore pari al 10% della deviazione standard
ampiezza_rumore_uffici_2022 = devstd_2022_uffici*0.1;
ampiezza_rumore_uffici_2021 = devstd_2021_uffici*0.1;
ampiezza_rumore_irraggiamento_2022 = devstd_2022_irraggiamento*0.1;
ampiezza_rumore_irraggiamento_2021 = devstd_2021_irraggiamento*0.1;

%% aggiunta del rumore
%rumore = intensita_rumore * randn(size(A, 1), 1);

forecast_uffici_2022 = data_2022(:,1) + ampiezza_rumore_uffici_2022*randn(size(data_2022,1),1);

for i=(1:size(data_2022,1))
    if data_2022(i,4)~=0
        forecast_irraggiamento_2022(i) = data_2022(i,4) + ampiezza_rumore_irraggiamento_2022*randn(1,1);
    else
        forecast_irraggiamento_2022(i) = 0;
    end
end
forecast_irraggiamento_2022 = forecast_irraggiamento_2022'.*(forecast_irraggiamento_2022'>=0);

%% Creazione dei file da passare al codice
ore = [1:24]';
ore_anno = repmat(ore, 365, 1);

Ir = [ore_anno,forecast_irraggiamento_2022,data_2022(:,4)];
Uffici = [ore_anno,forecast_uffici_2022,data_2022(:,1)];

