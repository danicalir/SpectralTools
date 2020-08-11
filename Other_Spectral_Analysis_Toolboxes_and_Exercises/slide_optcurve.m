function [kopt,err,optcurve] = slide_optcurve(k,mapped,cls,cns,ils,ins)

%[kopt,err,optcurve] = slide_optcurve(k,mapped,cls,cns,ils,ins);
%Computes the error index between two landslide maps following Carrara et
%al. (1992).  
%
%Inputs from slide_stats.m.  
%
%kopt = optimal cutoff spectral power sum and counts of correctly and
%incorrectly identified landslide and non-landslide nodes
%err = minimum error index
%optcurve = vector of error indices for cutoff values in k
%
%A.M. Booth (updated 11/2008)

%Count number of landslide nodes:
als = nansum(nansum(mapped,2),1);
als = als*ones(1,length(k));

%Compute the error index and find minimum:
optcurve = (1 - cls./(als + ils));
[err,ind] = min(optcurve);
kopt = [k(ind) cls(ind) cns(ind) ils(ind) ins(ind)];