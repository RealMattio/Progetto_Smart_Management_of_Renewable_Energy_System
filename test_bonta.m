% Dati di esempio
Ir = xlsread('Hybrid_model_single_forecast_irragiamento96h.xlsx');
% Add zeros rows and exchange the columns for idexes compatibility
Ir = [zeros(4440,3); Ir(:,1), Ir(:,3), Ir(:,2); zeros(2976,3)];
% Test della bontà delle stime
%Ir = [Ir(:,1), Ir(:,2), Ir(:,2)];

%load FALSA_previsione_irraggiamento.mat
% Format:[hour,forecasted Ir (°C), actual Ir (°C)] - Variable name: Ir



% Uncontrolled loads
% Format: [step, real power [W], forecasted power [W]]
Uffici = xlsread("Hybrid_model_single_forecast_24h.xlsx");
% Add zeros rows and exchange the columns for idexes compatibility
Uffici = [zeros(4368,3); Uffici(:,1), Uffici(:,3), Uffici(:,2); zeros(2952,3)];
% Test della bontà delle stime
%Uffici = [Uffici(:,1), Uffici(:,2), Uffici(:,2)];
%load FALSA_previsione_uffici.mat
% Format: [hour, forecasted power [W], actual power [W]] - Variable name : Uffici




predictions = Uffici(:,2);
real_data = Uffici(:,3);

% Calcola metriche e genera dati rumorosi
[MAE, RMSE, real_data_noisy] = calculate_metrics(predictions, real_data);

Uffici_worst = [Uffici(:,1), Uffici(:,3), real_data_noisy];
writematrix(Uffici_worst, 'Previsione_uffici_peggiore.xlsx');



function [MAE, RMSE, real_data_noisy] = calculate_metrics(predictions, real_data)
    % Verifica che gli input siano vettori della stessa lunghezza
    if ~isvector(predictions) || ~isvector(real_data)
        error('Gli input devono essere vettori.');
    end
    if length(predictions) ~= length(real_data)
        error('I vettori devono avere la stessa lunghezza.');
    end
    
    % Converti i vettori in colonne
    predictions = predictions(:);
    real_data = real_data(:);
    
    % Calcola MAE (Mean Absolute Error)
    MAE = mean(abs(predictions - real_data));
    
    % Calcola RMSE (Root Mean Squared Error)
    RMSE = sqrt(mean((predictions - real_data).^2));
    
    % Genera il vettore con rumore per peggiorare i dati reali
    error = predictions - real_data;

    valore = randi([0, 1], length(error),1) * 2 - 1;
    real_data_noisy = real_data - valore.*error;  % Aggiunge rumore per raddoppiare l'errore
    
    % Opzionale: Verifica che i nuovi errori siano effettivamente raddoppiati
    % MAE_new = mean(abs(predictions - real_data_noisy));
    % RMSE_new = sqrt(mean((predictions - real_data_noisy).^2));
    % fprintf('MAE originale: %.4f, MAE nuovo: %.4f\n', MAE, MAE_new);
    % fprintf('RMSE originale: %.4f, RMSE nuovo: %.4f\n', RMSE, RMSE_new);
end