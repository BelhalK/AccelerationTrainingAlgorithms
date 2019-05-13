function f = bolus_3cpt_Vkk12k21k13k31( phi, Id, X )
%Regression function for a PK model.
%
%Model information.
%         category: PK
%             name: bolus_3cpt_Vkk12k21k13k31
%      description:  IV bolus with three compartments and parameters: V, k, k12, k21, k13, k31
%                   %
% # regression var:  2
% # random effects:  6
% log_struct("L","N"): LLLL
%
%
% Function information. 
%       three input parameters 
%         phi (nI x 6) -  random effects. 
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


%% Initializing the random effect parameters.

V = phi(:,1); 
k = phi(:,2); 
k12 = phi(:,3); 
k21 = phi(:,4);
k13 = phi(:,5); 
k31 = phi(:,6);

a0 = k.*k21.*k31 ;
a1 = k.*k31 + k21.*k31 + k21.*k13 + k.*k21 + k31.*k12 ;
a2 = k + k12 + k13 + k21 + k31 ;
p = a1 - (a2.^2)./3 ;
q = (2.*(a2.^3)./27) - (a1.*a2./3) + a0 ;
r1 = sqrt(-(p.^3)./27) ;
r2 = 2.*(r1.^(1/3)) ;
p_phi = acos(-q./(2.*r1)) ./ 3 ;

alpha = -(cos(p_phi).*r2 - a2./3) ;
beta  = -(cos(p_phi+(2*pi/3)).*r2 - a2./3) ;
gamma = -(cos(p_phi+(4*pi/3)).*r2 - a2./3) ;

A = (1./V) .* ((k21-alpha)./(alpha-beta)) .* ((k31-alpha)./(alpha-gamma)) ;
B = (1./V) .* ((k21-beta)./(beta-alpha)) .* ((k31-beta)./(beta-gamma)) ;
C = (1./V) .* ((k21-gamma)./(gamma-beta)) .* ((k31-gamma)./(gamma-alpha)) ;

phi=[alpha, beta, gamma, A, B, C];

f = bolus_3cpt_alphabetagammaABC( phi, Id, X );