function f = bolus_1cpt_VCl( phi, Id, X )
%Regression function for a PK model.
%
%Model information.
%         category: PK
%             name: bolus_1cpt_VCl
%      description:  IV bolus with one compartment and parameters: V, Cl
%                   %
% # regression var:  2
% # random effects:  2
% log_struct("L","N"): LL
%
%
% Function information. 
%       three input parameters 
%         phi (nI x 2) -  random effects. 
%           All the random effects are Log-Normal.
%
%         X   struct with regression data. Fields
%              type_dose  - String with the type of the doses (single dose "sd", multiple dose "md", steady state "ss"). 
%              obs:    (nO x 1)- regression variables. 
%              dose:   (nD x ?) - doses informations (t_dose, dose, ...). 
%              t_inf:  (nO x 1) - Number of doses before the observation (only for "md"). 
%
%         Id   struct with individuals index informations. 
%              obs:    (nO x 1) - individual index for each observation. 
%             dose:    (nD x 1) - individual index for each dose. 
%              Nbs:    (nI x 2) - Number of doses and observations per individual. 
%
% where 
%     nI -> Number of individuals.
%     nO -> Number of observations.
%     nD -> Number of doses informations ( nD=nI for "sd" and "ss" )
%
%

phi(:,2)=phi(:,2)./phi(:,1);  % k = Cl/V
f = bolus_1cpt_Vk( phi, Id, X );