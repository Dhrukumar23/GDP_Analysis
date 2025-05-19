%% import the data

france = readtable("France GDP 1995 -2024.csv");
Germany = readtable("Germany  GDP 1995 2024.csv");
Italy = readtable ("Italy  GDP 1995 - 2024.csv");
Portugal = readtable("Portugal  GDP 1995 -2024.csv");
Spain =  readtable ("Spain  GDP 1995 -2024.csv");

data_f = table2array(france(:,2));
data_G = table2array(Germany(:,2));
data_I = table2array(Italy(:,2));
data_P = table2array(Portugal(:,2));
data_S = table2array(Spain(:,2));
year = table2array(Spain(:,1));

% PLOT OF THE TIME SERIES DATA 
figure;
plot(year,data_f);
title('Constanat GDP per Capita of 5 countires from 01/01/1995 to 01/01/2024');
xlabel('YEARS');
ylabel('USD$');

hold on
plot(year,data_G);
plot(year,data_I);
plot(year,data_P);
plot(year,data_S);
legend('france','Germany','Italy','Portugal','Spain');
colorbar;

% GROTH RATA
datal_f = diff(data_f)./data_f(1:end-1)*100;
datal_G = diff(data_G)./data_G(1:end-1)*100;
datal_I= diff(data_I)./data_I(1:end-1)*100;
datal_P = diff(data_P)./data_P(1:end-1)*100;
datal_S = diff(data_S)./data_S(1:end-1)*100;
ptime = year(2:end);

figure;
plot(ptime,datal_f,'-*');
title('The Constant GDP Rate of change per capita')
xlabel('YEARS')
ylabel('Rate of change')

hold on 
plot (ptime,datal_G,'-');
plot (ptime,datal_I,'-p');
plot (ptime,datal_P,'-d');
plot (ptime,datal_S,'-*');
legend('france','Germany','Italy','Portugal','Spain')
colorbar;

% ploting of histograph for each country
figure;
tiledlayout('flow');
nexttile
nbins = 100;
histogram(datal_f,nbins);
title('Histogaraph of farnce GDP')
xlabel ('% rate of change');
nexttile
nbins=100;
histogram(datal_G,nbins);
title('Histogaraph of Germany GDP')
xlabel ('% rate of change');
nexttile
nbins=100;
histogram(datal_I,nbins);
title('Histogaraph of Italy GDP')
xlabel ('% rate of change');
nexttile
nbins=100;
histogram(datal_P,nbins);
title('Histogaraph of Portugel GDP')
xlabel ('% rate of change');
nexttile
nbins=100;
histogram(data_S,nbins);
title('Histogaraph of Spain GDP')
xlabel ('% rate of change');

% loction estiomater

% computattion of mean for the GDP annual rata of change 
fprintf('\n');
france_mean = mean(datal_f);
Germany_mean = mean(datal_G);
Italy_mean = mean(datal_I);
Portugal_mean = mean(datal_P);
Spain_mean = mean(datal_S);

fprintf( 'the mean of France GDP rate of change is : %f\n',france_mean );
fprintf( 'the mean of Germany GDP rate of change is : %f\n',Germany_mean );
fprintf( 'the mean of Italy GDP rate of change is : %f\n',Italy_mean );
fprintf( 'the mean of Portugal GDP rate of change is : %f\n',Portugal_mean );
fprintf( 'the mean of Spain GDP rate of change is : %f\n ',Spain_mean );

%computattion of median for thr GDP annual rate of change in percentage 

fprintf('\n');
france_median = median(datal_f);
Germany_median = median(datal_G);
Italy_median = median(datal_I);
Portugal_median = median(datal_P);
Spain_median = median(datal_S);

fprintf( 'the median of France GDP rate of change is : %f\n',france_median );
fprintf( 'the median of Germany GDP rate of change is : %f\n',Germany_median );
fprintf( 'the median of Italy GDP rate of change is : %f\n',Italy_median );
fprintf( 'the median of Portugal GDP rate of change is : %f\n',Portugal_median );
fprintf( 'the median of Spain GDP rate of change is : %f\n',Spain_median );

%computattion of modefor thr GDP annual rate of change in percentage 

france_mode = mode(datal_f);
Germany_mode = mode(datal_G);
Italy_mode = mode(datal_I);
Portugal_mode = mode(datal_P);
Spain_mode = mode(datal_S);

fprintf( 'the mode of France GDP rate of change is : %f\n',france_mode );
fprintf( 'the mode of Germany GDP rate of change is : %f\n',Germany_mode );
fprintf( 'the mode of Italy GDP rate of change is : %f\n',Italy_mode );
fprintf( 'the mode of Portugal GDP rate of change is : %f\n',Portugal_mode );
fprintf( 'the mode of Spain GDP rate of change is : %f\n',Spain_mode );

tablea=table({'france';'germany';'italy';'portugal';'spain'},...
    [france_mean;Germany_mean;Italy_mean;Portugal_mean;Spain_mean],...
    [france_median;Germany_median;Italy_median;Portugal_median;Spain_mean],...
    [france_mode;Germany_mode;Italy_mode;Portugal_mode;Spain_mode],...
    'VariableName',{'country','mean','median','mode'});
disp(tablea);

%Disperion estimator (compution of GDP rate of change ranges)
fprintf('\n');
france_range = max(datal_f)-min(datal_f);
Germany_range = max(datal_G)-min(datal_G);
Italy_range = max(datal_I)-min(datal_I);
Portugal_range = max(datal_P)-min(datal_P);
Spain_range = max(datal_S)-min(datal_S);

fprintf( 'the range of France GDP rate of change is : %f\n',france_range);
fprintf( 'the range of Germany GDP rate of change is : %f\n',Germany_range );
fprintf( 'the range of Italy GDP rate of change is : %f\n',Italy_range );
fprintf( 'the range of Portugal GDP rate of change is : %f\n',Portugal_range );
fprintf( 'the range of Spain GDP rate of change is : %f\n',Spain_range );

%computation of constant GDP interquartiles 
fprintf('\n');
france_iqr = quantile(datal_f,0.75)-quantile(datal_f,0.25);
Germany_iqr = quantile(datal_f,0.75)-quantile(datal_G,0.25);
Italy_iqr = quantile(datal_f,0.75)-quantile(datal_I,0.25);
Portugal_iqr = quantile(datal_f,0.75)-quantile(datal_P,0.25);
Spain_iqr = quantile(datal_f,0.75)-quantile(datal_S,0.25);

fprintf( 'the interquartile of France GDP rate of change is : %f\n',france_iqr);
fprintf( 'the interquartile of Germany GDP rate of change is : %f\n',Germany_iqr );
fprintf( 'the interquartile of Italy GDP rate of change is : %f\n',Italy_iqr );
fprintf( 'the interquartile of Portugal GDP rate of change is : %f\n',Portugal_iqr );
fprintf( 'the interquartile of Spain GDP rate of change is : %f\n',Spain_iqr );

% Box polt %
figure;
boxplot([datal_f,datal_G,datal_I,datal_P,datal_S],'Labels',{'France','Garmany','Italy','Portugal','Spain'})
ylabel('% Gdp rate of change')

% Standard deviation 
fprintf('\n');
france_std = std(datal_f);
Germany_std = std(datal_G);
Italy_std = std(datal_I);
Portugal_std = std(datal_P);
Spain_std = std(datal_S);

fprintf( 'the Standard deviation of France GDP rate of change is : %f\n',france_std);
fprintf( 'the Standard deviation  of Germany GDP rate of change is : %f\n',Germany_std);
fprintf( 'the Standard deviation  of Italy GDP rate of change is : %f\n',Italy_std);
fprintf( 'the Standard deviation  of Portugal GDP rate of change is : %f\n',Portugal_std );
fprintf( 'the Standard deviation  of Spain GDP rate of change is : %f\n',Spain_std );

%Variance 
fprintf('\n');
France_var = var(datal_f);
Garmany_var = var(datal_G);
Italy_var = var(datal_I);
Portugal_var = var(datal_P);
Spain_var = var(datal_S);

fprintf( 'the Sample Variance of France GDP rate of change is : %f\n',France_var);
fprintf( 'the   Sample Variance of Germany GDP rate of change is : %f\n',Garmany_var );
fprintf( 'the  Sample Variance  of Italy GDP rate of change is : %f\n',Italy_var);
fprintf( 'the  Sample Variance  of Portugal GDP rate of change is : %f\n',Portugal_var );
fprintf( 'the  Sample Variance  of Spain GDP rate of change is : %f\n',Spain_var);

D_var_std=table({'france';'germany';'italy';'portugal';'spain'},...
    [France_var;Garmany_var;Italy_var;Portugal_var;Spain_var],...
    [france_std;Germany_std;Italy_std;Portugal_std;Spain_std],...
    'VariableName',{'country','variance','Standard deviation'});
disp(D_var_std);


% std normal disttibuton 
figure;
std_rates_France = (datal_f - france_mean)./std(datal_f);
std_rates_France;
std_rates_Germany = (datal_G- Germany_mean)./std(datal_G);
std_rates_Germany;
std_rates_Italy = (datal_I-Italy_mean)./std(datal_I);
std_rates_Italy;
std_rates_Portugal = (datal_P - Portugal_mean)./std(datal_P);
std_rates_Portugal;
std_rates_Spain = (datal_S- Spain_mean)./std(datal_S);
std_rates_Spain;

% Histogram and comparison for Each Countries
sgtitle('Comparison of Standardized Rates with Standard Normal Distribution');
No_Bins = 20;
   % Histogram and comparison for france
   % Subplot for France (Normalize to PDF)
subplot(3, 2, 1);
histogram(std_rates_France, No_Bins, 'Normalization', 'pdf', 'FaceColor', 'y'); 
hold on;
x = linspace(-4, 4, 100);
y = normpdf(x, 0, 1); 
plot(x, y, 'r', 'LineWidth', 3); 
title('Histogram and comparison for france');
xlabel('Standardized Rates');
ylabel('Probability Density');
legend('standardized Growth rate', 'Standard Normal','Location', 'Best');
grid on;

sgtitle('Comparison of Standardized Rates with Standard Normal Distribution');
No_Bins = 20;
    % Histogram and comparison for Germany
        % Subplot for Germany (Normalize to PDF)
subplot(3, 2, 2);
histogram(std_rates_Germany, No_Bins, 'Normalization', 'pdf', 'FaceColor', 'b'); 
hold on;
x = linspace(-4, 4, 100);
y = normpdf(x, 0, 1); 
plot(x, y, 'r', 'LineWidth', 3); 
title('Histogram and comparison for Germany');
xlabel('Standardized Rates');
ylabel('Probability Density');
legend('standardized Growth rate', 'Standard Normal','Location', 'Best');
grid on;

sgtitle('Comparison of Standardized Rates with Standard Normal Distribution');
No_Bins = 20;
    % Histogram and comparison for Italy
        % Subplot for Italy (Normalize to PDF)
subplot(3, 2, 3);
histogram(std_rates_Italy, No_Bins, 'Normalization', 'pdf', 'FaceColor', 'b'); 
hold on;
x = linspace(-4, 4, 100);
y = normpdf(x, 0, 1); 
plot(x, y, 'r', 'LineWidth', 3); 
title('Histogram and comparison for Italy');
xlabel('Standardized Rates');
ylabel('Probability Density');
legend('standardized Growth rate', 'Standard Normal','Location', 'Best');
grid on;

sgtitle('Comparison of Standardized Rates with Standard Normal Distribution');
No_Bins = 20;
    % Histogram and comparison for Portugal
        % Subplot for Portugal (Normalize to PDF)
subplot(3, 2, 4);
histogram(std_rates_Portugal, No_Bins, 'Normalization', 'pdf', 'FaceColor', 'y'); 
hold on;
x = linspace(-4, 4, 100);
y = normpdf(x, 0, 1); 
plot(x, y, 'r', 'LineWidth', 3); 
title('Histogram and comparison for Portugal');
xlabel('Standardized Rates');
ylabel('Probability Density');
legend('standardized Growth rate', 'Standard Normal','Location', 'Best');
grid on;

sgtitle('Comparison of Standardized Rates with Standard Normal Distribution');
No_Bins = 20;
    % Histogram and comparison for Spain
        % Subplot for Spain (Normalize to PDF)
subplot(3, 2, 5);
histogram(std_rates_Spain, No_Bins, 'Normalization', 'pdf', 'FaceColor', 'y'); 
hold on;
x = linspace(-4, 4, 100);
y = normpdf(x, 0, 1); 
plot(x, y, 'r', 'LineWidth', 3); 
title('Histogram and comparison for Spain');
xlabel('Standardized Rates');
ylabel('Probability Density');
legend('standardized Growth rate', 'Standard Normal','Location', 'Best');
grid on;

% Dependence estimators 
% Scatter plot %
figure;
tiledlayout("flow");
nexttile;
plot(datal_f,datal_G,'*');
title('scatter plot 1')
xlabel('france');
ylabel('Germany')
nexttile;
plot(datal_f,datal_I,'*');
title('scatter plot 2')
xlabel('france');
ylabel('Italy');
nexttile;
plot(datal_f,datal_P,'*');
title('scatter plot 3')
xlabel('france');
ylabel('Portugal');
nexttile;
plot(datal_f,datal_S,'*');
title('scatter plot 4')
xlabel('france');
ylabel('Spain');
nexttile;
plot(datal_G,datal_I,'*');
title('scatter plot 5')
xlabel('Germany');
ylabel('Italy');
nexttile;
plot(datal_G,datal_P,'*');
title('scatter plot 6')
xlabel('Germany');
ylabel('Portugal');
nexttile;
plot(datal_G,datal_S,'*');
title('scatter plot 7')
xlabel('Germany');
ylabel('Spain');
nexttile;
plot(datal_I,datal_P,'*');
title('scatter plot 8')
xlabel('Italy');
ylabel('Portugal');
nexttile;
plot(datal_I,datal_S,'*');
title('scatter plot 9')
xlabel('Italy');
ylabel('Spain');
nexttile;
plot(datal_P,datal_S,'*');
title('scatter plot 10')
xlabel('Portugal');
ylabel('Spain');

%correlation
cor_france_Germany = corr([datal_f,datal_G]);
cor_france_Italy = corr([datal_f,datal_I]);
cor_france_Portugal = corr([datal_f,datal_P]);
cor_france_Spain = corr([datal_f,datal_S]);
cor_Germany_Itlay = corr([datal_G,datal_I]);
cor_Germany_Portugal = corr([datal_G,datal_P]);
cor_Germany_Spain = corr([datal_G,datal_S]);
cor_Italy_Portugal = corr([datal_I,datal_P]);
cor_Italy_Spain = corr([datal_I,datal_S]);
cor_Portugal_Spain = corr([datal_P,datal_S]);

fprintf('\n')
fprintf('cor france is Germany:%f\n',cor_france_Germany(2));
fprintf('cor france is Italy:%f\n',cor_france_Italy(2));
fprintf('cor france is Portugal:%f\n',cor_france_Portugal(2));
fprintf('cor france is Spain:%f\n',cor_france_Spain(2));
fprintf('cor Germany is Italy:%f\n',cor_Germany_Itlay(2));
fprintf('cor Germany is Portugal:%f\n',cor_Germany_Portugal(2));
fprintf('cor Germany is Spain:%f\n',cor_Germany_Spain(2));
fprintf('cor Italy is Portugal:%f\n',cor_Italy_Portugal(2));
fprintf('cor Italy is Spain:%f\n',cor_Italy_Spain(2));
fprintf('cor Portugel is Spain:%f\n',cor_Portugal_Spain(2));


cov_france_Germany = cov([datal_f,datal_G]);
cov_france_Italy = cov([datal_f,datal_I]);
cov_france_Portugal = cov([datal_f,datal_P]);
cov_france_Spain = cov([datal_f,datal_S]);
cov_Germany_Itlay = cov([datal_G,datal_I]);
cov_Germany_Portugal = cov([datal_G,datal_P]);
cov_Germany_Spain = cov([datal_G,datal_S]);
cov_Italy_Portugal = cov([datal_I,datal_P]);
cov_Italy_Spain = cov([datal_I,datal_S]);
cov_Portugal_Spain = cov([datal_P,datal_S]);
 
fprintf('\n')
fprintf('cov france is Germany:%f\n',cov_france_Germany(2))
fprintf('cov france is Italy:%f\n',cov_france_Italy(2));
fprintf('cov france is Portugal:%f\n',cov_france_Portugal(2));
fprintf('cov france is Spain:%f\n',cov_france_Spain(2));
fprintf('cov Germany is Italy:%f\n',cov_Germany_Itlay(2));
fprintf('cov Germany is Portugal:%f\n',cov_Germany_Portugal(2));
fprintf('cov Germany is Spain:%f\n',cov_Germany_Spain(2));
fprintf('cov italy is Portugal:%f\n',cov_Italy_Portugal(2));
fprintf('cov Italy is Spain:%f\n',cov_Italy_Spain(2));
fprintf('cov Portugel is Spain:%f\n',cov_Portugal_Spain(2));

T_cor_cov= table({'france_Germany';'france_Italy';'france_Portugal';'france_Spain';'Germany_Itlay';'Germany_Portugal';'Germany_Spain';'Italy_Portugal';'Italy_Spain';'Portugal_Spain'},...
    [cor_france_Germany(2);cor_france_Italy(2);cor_france_Portugal(2);cor_france_Spain(2);cor_Germany_Itlay(2);cor_Germany_Portugal(2);cor_Germany_Spain(2);cor_Italy_Portugal(2);cor_Italy_Spain(2);cor_Portugal_Spain(2)],...
    [cov_france_Germany(2);cov_france_Italy(2);cov_france_Portugal(2);cov_france_Spain(2);cov_Germany_Itlay(2);cov_Germany_Portugal(2);cov_Germany_Spain(2);cov_Italy_Portugal(2);cov_Italy_Spain(2);cov_Portugal_Spain(2)],....
    'VariableName',{'country','COR','COV'});
disp(T_cor_cov);

%COMPUTATION OF CONFIDENCE INTERVALS FOR EACH COUNTRY
% I STARTED BY DETEMINING THE NUMBER OF OBSERVATION (NOBS) THEN USED
% STUDENT-T TO COMPUTE CONFIDENCE INTERVALS

nobs = size (datal_f,1);

%Confidence interval France 
ci_f90 = [france_mean-tinv(.95,nobs-1)*france_std/sqrt(nobs),france_mean+tinv(.95,nobs-1)*france_std/sqrt(nobs)];
ci_f95 = [france_mean-tinv(.975,nobs-1)*france_std/sqrt(nobs),france_mean+tinv(.975,nobs-1)*france_std/sqrt(nobs)];
ci_f99 = [france_mean-tinv(.995,nobs-1)*france_std/sqrt(nobs),france_mean+tinv(.995,nobs-1)*france_std/sqrt(nobs)];

 fprintf('\n')
 fprintf('the for france ci90is%f\n',ci_f90);
 fprintf('the for france ci95is%f\n',ci_f95);
 fprintf('the for france ci99is%f\n',ci_f99);

%Confidence interval Germany
ci_g90 = [Germany_mean-tinv(.95,nobs-1)*Germany_std/sqrt(nobs),Germany_mean+tinv(.95,nobs-1)*Germany_std/sqrt(nobs)];
ci_g95 = [Germany_mean-tinv(.975,nobs-1)*Germany_std/sqrt(nobs),Germany_mean+tinv(.975,nobs-1)*Germany_std/sqrt(nobs)];
ci_g99 = [Germany_mean-tinv(.995,nobs-1)*Germany_std/sqrt(nobs),Germany_mean+tinv(.997,nobs-1)*Germany_std/sqrt(nobs)];
 
fprintf('\n')
 fprintf('the for Germany ci90is%f\n',ci_g90);
 fprintf('the for  Germany ci95is%f\n',ci_g95);
 fprintf('the for  Germany ci99is%f\n',ci_g99);

%Confidence interval Italy 
ci_I90 = [Italy_mean-tinv(.95,nobs-1)*Italy_std/sqrt(nobs),Italy_mean+tinv(.95,nobs-1)*Italy_std/sqrt(nobs)];
ci_I95 = [Italy_mean-tinv(.975,nobs-1)*Italy_std/sqrt(nobs),Italy_mean+tinv(.975,nobs-1)*Italy_std/sqrt(nobs)];
ci_I99 = [Italy_mean-tinv(.995,nobs-1)*Italy_std/sqrt(nobs),Italy_mean+tinv(.995,nobs-1)*Italy_std/sqrt(nobs)];

fprintf('\n')
 fprintf('the for Italy ci90is%f\n',ci_I90);
 fprintf('the for  Italy ci95is%f\n',ci_I95);
 fprintf('the for  Italy ci99is%f\n',ci_I99);

%Confidence interval Portugal
ci_P90 = [Portugal_mean-tinv(.95,nobs-1)*Portugal_std/sqrt(nobs),Portugal_mean+tinv(.95,nobs-1)*Portugal_std/sqrt(nobs)];
ci_P95 = [Portugal_mean-tinv(.975,nobs-1)*Portugal_std/sqrt(nobs),Portugal_mean+tinv(.975,nobs-1)*Portugal_std/sqrt(nobs)];
ci_P99 = [Portugal_mean-tinv(.995,nobs-1)*Portugal_std/sqrt(nobs),Portugal_mean+tinv(.995,nobs-1)*Portugal_std/sqrt(nobs)];

fprintf('\n')
fprintf('the for Portugal ci90is%f\n',ci_P90);
 fprintf('the for  Portugal ci95is%f\n',ci_P95);
fprintf('the for  Portugal ci99is%f\n',ci_P99);

%Confidence interval Spain
ci_S90 = [Spain_mean-tinv(.95,nobs-1)*Spain_std/sqrt(nobs),Spain_mean+tinv(.95,nobs-1)*Spain_std/sqrt(nobs)];
ci_S95 = [Spain_mean-tinv(.975,nobs-1)*Spain_std/sqrt(nobs),Spain_mean+tinv(.975,nobs-1)*Spain_std/sqrt(nobs)];
ci_S99 = [Spain_mean-tinv(.995,nobs-1)*Spain_std/sqrt(nobs),Spain_mean+tinv(.995,nobs-1)*Spain_std/sqrt(nobs)];

fprintf ('\n')
fprintf('the for Spain ci90is%f\n',ci_S90);
fprintf('the for  Spain ci95is%f\n',ci_S95);
fprintf('the for  Spain ci99is%f\n',ci_S99);

