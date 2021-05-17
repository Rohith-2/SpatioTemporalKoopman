
clear all
clc
% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Open", "Var3", "Var4", "Close", "Var6", "Volume"];
opts.SelectedVariableNames = ["Open", "Close", "Volume"];
opts.VariableTypes = ["string", "double", "string", "string", "double", "string", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var3", "Var4", "Var6"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var3", "Var4", "Var6"], "EmptyFieldRule", "auto");

% Import the data
BTCUSD = readtable("/Volumes/Rohith/SpatioTemporalKoopman/data/BTC-USD.csv", opts);
DOGEUSD = readtable("/Volumes/Rohith/SpatioTemporalKoopman/data/DOGE-USD.csv", opts);
ETHUSD = readtable("/Volumes/Rohith/SpatioTemporalKoopman/data/ETH-USD.csv", opts);

% Convert to output type
BTC = table2array(BTCUSD)';
DOGE = table2array(DOGEUSD)';
ETH = table2array(ETHUSD)';

% Clear temporary variables
clear opts
data(:,:,1) = (ETH);
data(:,:,2) = (DOGE);
data(:,:,3) = (BTC);

%Removal of NaN's with previous days trade
for k = 1:3
    for i = 1:2111
        for j=1:3
            if(isnan(data(j,i,k)))
                data(j,i,k) = data(j,i-1,k);
            end

        end
    end
end