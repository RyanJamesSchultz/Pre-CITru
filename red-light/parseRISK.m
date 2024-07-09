function [Rd,Rc,Ri,Rf,Re,M]=parseRISK(filename)
  % Simple routine to parse csv files from OQ.
  
  % Load in the data.
  % Canton, displaced, content loss [CHF], injuries, fatalities, economic loss [CHF], % displaced, mag
  command=['cat ', filename, ' | ', ... % Load in all the EQ data files.
           'awk -F, ''/Canton/{next} {print $2,$3,$4,$5,$6,$7,$8}'' > temp.RiskParse ']; % Parse csv files, and put in NaN for null values.
  system(command);
  data=load('temp.RiskParse');
  system('rm -f temp.RiskParse');
  Rd=data(:,1); % Displaced.
  Rc=data(:,2); % Content loss [CHF].
  Ri=data(:,3); % Injuries.
  Rf=data(:,4); % Fatalities.
  Re=data(:,5); % Economic loss (structural & non-structural) [CHF].
  M=data(:,7);
  
return


