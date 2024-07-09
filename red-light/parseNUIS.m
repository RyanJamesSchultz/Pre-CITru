function [N1,N2,N3,N4,N5,N6,M]=parseNUIS(filename)
  % Simple routine to parse csv files from OQ.
  
  % Load in the data.
  % Canton, displaced, content loss [CHF], injuries, fatalities, economic loss [CHF], % displaced, mag
  command=['cat ', filename, ' | ', ... % Load in all the EQ data files.
           'awk -F, ''/Canton/{next} {print $2,$3,$4,$5,$6,$7,$8}'' > temp.RiskParse ']; % Parse csv files, and put in NaN for null values.
  system(command);
  data=load('temp.RiskParse');
  system('rm -f temp.RiskParse');
  N1=data(:,1); % CDI1.
  N2=data(:,2); % CDI2.
  N3=data(:,3); % CDI3.
  N4=data(:,4); % CDI4.
  N5=data(:,5); % CDI5.
  N6=data(:,6); % CDI6.
  M=data(:,7);
  
return


