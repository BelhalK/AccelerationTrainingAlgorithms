function f = bolus_3cpt_alphabetagammaABC( phi, Id, X )
%Regression function for a PK model.
%
%Model information.
%         category: PK
%             name: bolus_3cpt_alphabetagammaABC
%      description:  IV bolus with three compartments and parameters: Alpha, Beta, Gamma, A, B, C
%
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

alpha=phi(Id.obs,1);
beta=phi(Id.obs,2);
gamma=phi(Id.obs,3);
A=phi(Id.obs,4);
B=phi(Id.obs,5);
C=phi(Id.obs,6);

TAD=X.TAD;%Time after the last dose

%Here you can evaluate your model function "f" for each type of model "sd","md","ss"
%

switch X.type_dose
    case'sd'
        %     Single dose model
        dt    = TAD; % difference between time of dose and time of observation
        D     = X.dose(Id.obs, 2);  %  dose
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Evaluate here your function for single dose
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f     = D.*(A.*exp(-alpha.*dt)+B.*exp(-beta.*dt)+C.*exp(-gamma.*dt));
        f     = f.*(dt>=0);
        %
        
    case 'md'
        %     Multiple dose model
        D        = X.dose(:,2);  % doses
        %
        alpha_dose=alpha(Id.obs_dose);
        beta_dose=beta(Id.obs_dose);
        gamma_dose=gamma(Id.obs_dose);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Evaluate here your functions (for multiple doses).
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ATAD=A.*exp(-alpha.*(TAD));
        BTAD=B.*exp(- beta.*(TAD));
        CTAD=C.*exp(-gamma.*(TAD));
        
        %time between doses
        TBD=X.TBD;
        
        ATBD=X.B_dose;
        BTBD=ATBD;CTBD=ATBD;
        
        ATBD(X.I_dose)=exp(-alpha_dose.*(TBD));
        BTBD(X.I_dose)=exp(- beta_dose.*(TBD));
        CTBD(X.I_dose)=exp(-gamma_dose.*(TBD));
        
        %    ii=(2:numel(t_dose))';
        B1=X.B_dose;
        B1(X.I_dose)=D;
        B2=B1;
        B3=B1;
        
        for d=3:size(B1,1),
            B1(d,:)=B1(d,:) + B1(d-1,:).*ATBD(d,:);
            B2(d,:)=B2(d,:) + B2(d-1,:).*BTBD(d,:);
            B3(d,:)=B3(d,:) + B3(d-1,:).*CTBD(d,:);
        end
        f=ATAD.*B1(X.J_dose)+BTAD.*B2(X.J_dose)+CTAD.*B3(X.J_dose);
        
        %%
    otherwise
        %     Steady state model
        dt    = TAD; % difference between time of dose and time of observation
        D     = X.dose(Id.obs,2);  %  dose
        tau   = X.dose(Id.obs,3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Evaluate here your function (for steady state)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        D1=max(eps,1-exp(-alpha.*tau));
        D2=max(eps,1-exp(-beta.*tau));
        D3=max(eps,1-exp(-gamma.*tau));
        
        f =D.*(A.*exp(-alpha.*dt)./D1+B.*exp(-beta.*dt)./D2+C.*exp(-gamma.*dt)./D3);
        
        f=f.*(dt>=0);
end
f=max(0,f);
