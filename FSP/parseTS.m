function [vID,vP,tIDs]=parseTS(filename)
  % Simple routine to parse csv files from OQ.
  
  % Load in the data.
  command1=['cat ', filename, ' | ', ... % Load in all the EQ data files.
           'awk ''/TRGL/{print $2,$3,$4}'' > temp.TStParse ']; % Parse csv files, and put in NaN for null values.
  system(command1);
  TRGLdata=load('temp.TStParse');
  command2=['cat ', filename, ' | ', ... % Load in all the EQ data files.
           'awk ''/VRTX/{print $2,$3,$4,$5}'' > temp.TSvParse ']; % Parse csv files, and put in NaN for null values.
  system(command2);
  VRTXdata=load('temp.TSvParse');
  system('rm -f temp.TSvParse');
  system('rm -f temp.TStParse');
  vID=VRTXdata(:,1);
  vP=[VRTXdata(:,2), VRTXdata(:,3), VRTXdata(:,4)];
  tIDs=[TRGLdata(:,1), TRGLdata(:,2), TRGLdata(:,3)];

  return
end