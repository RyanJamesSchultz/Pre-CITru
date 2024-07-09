function [lat,lon,PGV,dPGV]=parseSM(filename)
  % Simple routine to parse csv files from ShakeMap.
  
  % Load in the data.
  % Canton, displaced, content loss [CHF], injuries, fatalities, economic loss [CHF], % displaced, mag
  command=['cat ', filename, ' | ', ... % Load in all the EQ data files.
           'awk -F, ''/lat/{next} {print $1,$2,$3,$4}'' > temp.SMParse ']; % Parse csv files, and put in NaN for null values.
  system(command);
  data=load('temp.SMParse');
  system('rm -f temp.SMParse');
  lat=data(:,1); % Displaced.
  lon=data(:,2); % Content loss [CHF].
  PGV=data(:,3); % PGV (m/s).
  dPGV=exp(data(:,4))/100; % uncertainty (m/s).
  
return


