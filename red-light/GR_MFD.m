function [PDF,CDF]=GR_MFD(m,Mb,params,norm_flag)
  % Function that computes the Gutenberg-Richter magnitude-frequency 
  % distribution (GR-MFD).
  % 
  % Written by Ryan Schultz.
  
  % Get the upper/lower magnitude bounds.
  m1=min(Mb);
  m2=max(Mb);
  
  % Get the parameters.
  a=params(1);
  b=params(2);
  
  % Define some useful values.
  f=b*log(10);
  x=10.^(-b*(m-m1));
  n=1-10.^(-b*(m2-m1));
  
  % The PDF & CDF for a bounded GR-MFD
  PDF=(f*x)/n; PDF(m>m2)=NaN; PDF(m<m1)=NaN;
  CDF=(1-x)/n; CDF(m>m2)=NaN; CDF(m<m1)=NaN;
  
  % Change the normalization, if flagged to.
  if(strcmpi(norm_flag,'count'))
      PDF=PDF*n*10^a;
      CDF=CDF*n*10^a;
  end
  
end