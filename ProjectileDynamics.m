classdef ProjectileDynamics < handle
    % TODO: Add tempStorage to derivative functions. That way not recomputing a ton of variables on each force model or Jacobian call

    properties
        projectile
        planet

        includeParamSTM = false;
    end

    properties (SetAccess = private)
        computeStateJacobian
    end
    

    methods
        % Constructor ==============================================================================

        function self = ProjectileDynamics(projectile, planet)
            if nargin ~= 2
                error("Not enough input arguments. Requires projectile and planet.")
            end

            self.projectile = projectile;
            self.planet = planet;
            
            switch Settings.DEFAULT_STATE_JACOBIAN_METHOD
                case "analytic"
                    self.computeStateJacobian = @self.computeAnalyticStateJacobian;
                
                case "numeric"
                    self.computeStateJacobian = @self.computeNumericStateJacobian;
            end
        end

        
        % Force model methods ======================================================================

        function F = computeGravityForce(self, state)
            theta = state(5);

            sinTheta = sin(theta);
            cosTheta = cos(theta);

            m = self.projectile.params(self.projectile.mIdx);

            g = self.planet.computeGravity();

            F = zeros(3, 1);
            F(1) = -m * g * sinTheta;
            F(3) =  m * g * cosTheta;
        end


        function [F, M] = computeAeroForce(self, state)
            z = state(3);
            theta = state(5);
            psi = state(6);
            u = state(7);
            v = state(8);
            w = state(9);
            p = state(10);
            q = state(11);
            r = state(12);

            cosTheta = cos(theta);
            sinTheta = sin(theta);
            cosPsi = cos(psi);
            sinPsi = sin(psi);

            d = self.projectile.params(self.projectile.dIdx);
            S = self.projectile.params(self.projectile.SIdx);
            deltaFins = self.projectile.params(self.projectile.deltaFinsIdx);
            
            h = -z;
            [rho, a] = self.planet.computeAtmosphere(h);
            [vWindx, vWindy] = self.planet.computeWind(h);

            uWind =  cosTheta * cosPsi * vWindx + cosTheta * sinPsi * vWindy;
            vWind = -sinPsi            * vWindx + cosPsi            * vWindy;
            wWind =  sinTheta * cosPsi * vWindx + sinTheta * sinPsi * vWindy;

            uAtm = u - uWind;
            vAtm = v - vWind;
            wAtm = w - wWind;
            
            VAtm = (uAtm ^ 2 + vAtm ^ 2 + wAtm ^ 2) ^ 0.5;
            
            sinAlphaCosBeta = wAtm / VAtm;
            sinBeta         = vAtm / VAtm;
            sinAlphaTotal   = (vAtm ^ 2 + wAtm ^ 2) ^ 0.5 / VAtm;

            mach = VAtm / a;
            qAtm = 0.5 * rho * VAtm ^ 2;

            pNormalized = p * d / VAtm;
            qNormalized = q * d / VAtm;
            rNormalized = r * d / VAtm;

            [CX0, CX2, CY0, CZ0, CNalpha0, CNalpha2, CNpalpha0, CNpalpha2, ...
             Cl0, Clp, Cldelta, Cm0, Cn0, CMalpha0, CMalpha2, CMpalpha0, CMpalpha2, CMq] = self.projectile.computeAeroCoeffs(mach);

            CX       = CX0       + CX2       * sinAlphaTotal ^ 2;
            CNalpha  = CNalpha0  + CNalpha2  * sinAlphaTotal ^ 2;
            CNpalpha = CNpalpha0 + CNpalpha2 * sinAlphaTotal ^ 2;
            CMalpha  = CMalpha0  + CMalpha2  * sinAlphaTotal ^ 2;
            CMpalpha = CMpalpha0 + CMpalpha2 * sinAlphaTotal ^ 2;

            CY = CY0 - CNalpha * sinBeta          + pNormalized * CNpalpha * sinAlphaCosBeta;
            CZ = CZ0 - CNalpha * sinAlphaCosBeta  - pNormalized * CNpalpha * sinBeta;
            
            Cl = Cl0 + pNormalized * Clp + Cldelta * deltaFins;
            Cm = Cm0 + CMalpha * sinAlphaCosBeta  + pNormalized * CMpalpha * sinBeta         + qNormalized * CMq;
            Cn = Cn0 - CMalpha * sinBeta          + pNormalized * CMpalpha * sinAlphaCosBeta + rNormalized * CMq;

            F = zeros(3, 1);
            F(1) = qAtm * S * CX;
            F(2) = qAtm * S * CY;
            F(3) = qAtm * S * CZ;

            M = zeros(3, 1);
            M(1) = qAtm * S * d * Cl;
            M(2) = qAtm * S * d * Cm;
            M(3) = qAtm * S * d * Cn;
        end


        % function dF_dv = computeAeroForcePartials(self, state)
        %     z = state(3);
        %     vx = state(4);
        %     vy = state(5);
        %     vz = state(6);
        % 
        %     S = self.projectile.params(self.projectile.SIdx);
        % 
        %     h = -z;
        %     [rho, a] = self.planet.computeAtmosphere(h);
        %     [vWindx, vWindy] = self.planet.computeWind(h);
        % 
        %     vAtmx = vx - vWindx;
        %     vAtmy = vy - vWindy;
        %     vAtmz = vz;
        %     VAtm = (vAtmx ^ 2 + vAtmy ^ 2 + vAtmz ^ 2) ^ 0.5;
        % 
        %     mach = VAtm / a;
        % 
        %     CD = self.projectile.computeAeroCoeffs(mach);
        % 
        %     k = -(rho * S * CD / 2);
        % 
        %     dF_dv = zeros(3, 3);
        % 
        %     dF_dv(1, 1) = k * (vAtmx ^ 2 / VAtm + VAtm);
        %     dF_dv(2, 1) = k * (vAtmx * vAtmy / VAtm);
        %     dF_dv(3, 1) = k * (vAtmx * vAtmz / VAtm);
        % 
        %     dF_dv(1, 2) = dF_dv(2, 1);
        %     dF_dv(2, 2) = k * (vAtmy ^ 2 / VAtm + VAtm);
        %     dF_dv(3, 2) = k * (vAtmy * vAtmz / VAtm);
        % 
        %     dF_dv(1, 3) = dF_dv(3, 1);
        %     dF_dv(2, 3) = dF_dv(3, 2);
        %     dF_dv(3, 3) = k * (vAtmz ^ 2 / VAtm + VAtm);
        % end

        
        % Derivative methods =======================================================================

        function stateDeriv = computeStateDeriv(self, state)
            nStates = self.projectile.nStates;

            theta = state(5);
            psi = state(6);
            u = state(7);
            v = state(8);
            w = state(9);
            p = state(10);
            q = state(11);
            r = state(12);

            cosTheta = cos(theta);
            sinTheta = sin(theta);
            tanTheta = sinTheta / cosTheta;
            secTheta = 1 / cosTheta;
            cosPsi = cos(psi);
            sinPsi = sin(psi);
            
            % TODO: Compute acceleration function?
            m = self.projectile.params(self.projectile.mIdx);
            Ix = self.projectile.params(self.projectile.IxxIdx);
            Iz = self.projectile.params(self.projectile.IzzIdx);
            dIxIz = Ix / Iz;

            FGrav = self.computeGravityForce(state);
            [FAero, M] = self.computeAeroForce(state);

            F = FGrav + FAero;

            stateDeriv = zeros(nStates, 1);

            stateDeriv(1) =  cosTheta * cosPsi * u - sinPsi * v + sinTheta * cosPsi * w;
            stateDeriv(2) =  cosTheta * sinPsi * u + cosPsi * v + sinTheta * sinPsi * w;
            stateDeriv(3) = -sinTheta          * u              + cosTheta          * w;

            stateDeriv(4) = p    + tanTheta * r;
            stateDeriv(5) =   q;
            stateDeriv(6) =        secTheta * r;

            stateDeriv(7) = F(1) / m         +            r * v -            q * w;
            stateDeriv(8) = F(2) / m - r * u                    - tanTheta * r * w;
            stateDeriv(9) = F(3) / m + q * u + tanTheta * r * v;

            stateDeriv(10) = M(1) / Ix;
            stateDeriv(11) = M(2) / Iz - dIxIz * r * p                    - tanTheta * r * r;
            stateDeriv(12) = M(3) / Iz + dIxIz * q * p + tanTheta * r * q;
        end


        function augStateDeriv = computeAugStateDeriv(self, augState)
            nStates = self.projectile.nStates;

            iStateEnd = nStates;
            iStateSTMEnd = nStates + nStates ^ 2;

            state = augState(1:iStateEnd);
            stateSTM = augState((iStateEnd + 1):iStateSTMEnd);
            if self.includeParamSTM
                paramSTM = augState((iStateSTMEnd + 1):end);
            end

            stateDeriv = self.computeStateDeriv(state);
            if self.includeParamSTM
                [stateSTMDeriv, paramSTMDeriv] = self.computeSTMDerivs(state, stateSTM, paramSTM);
                augStateDeriv = [stateDeriv; stateSTMDeriv; paramSTMDeriv];
            else
                stateSTMDeriv = self.computeSTMDerivs(state, stateSTM);
                augStateDeriv = [stateDeriv; stateSTMDeriv];
            end
        end


        function [stateSTMDeriv, paramSTMDeriv] = computeSTMDerivs(self, state, stateSTM, paramSTM)
            nStates = self.projectile.nStates;

            stateSTM = reshape(stateSTM, [nStates, nStates]);

            stateA = self.computeStateJacobian(state);

            stateSTMDeriv = stateA * stateSTM;
            stateSTMDeriv = stateSTMDeriv(:);

            if self.includeParamSTM
                nEstimatedParams = self.projectile.nEstimatedParams + self.planet.nEstimatedParams;

                paramSTM = reshape(paramSTM, [nStates, nEstimatedParams]);

                paramA = self.computeParamJacobian(state);

                paramSTMDeriv = stateA * paramSTM + paramA;
                paramSTMDeriv = paramSTMDeriv(:);
            end
        end


        % Jacobian methods =========================================================================

        function A = computeAnalyticStateJacobian(self, state)
            error("Analytic state Jacobian not yet implemented.")

            % nStates = self.projectile.nStates;
            % 
            % m = self.projectile.params(self.projectile.mIdx);
            % 
            % % dFGrav_dv = self.computeGravityForcePartials();
            % % dFAero_dv = self.computeAeroForcePartials(state);
            % % dF_ddv = dFGrav_dv + dFAero_dv;
            % 
            % dF_dv = self.computeAeroForcePartials(state);
            % 
            % A = zeros(nStates);
            % 
            % A(1, 4) = 1;
            % A(4, 4) = dF_dv(1, 1) / m;
            % A(5, 4) = dF_dv(2, 1) / m;
            % A(6, 4) = dF_dv(3, 1) / m;
            % 
            % A(2, 5) = 1;
            % A(4, 5) = dF_dv(1, 2) / m;
            % A(5, 5) = dF_dv(2, 2) / m;
            % A(6, 5) = dF_dv(3, 2) / m;
            % 
            % A(3, 6) = 1;
            % A(4, 6) = dF_dv(1, 3) / m;
            % A(5, 6) = dF_dv(2, 3) / m;
            % A(6, 6) = dF_dv(3, 3) / m;
        end


        function A = computeNumericStateJacobian(self, state)
            nStates = self.projectile.nStates;

            A = zeros(nStates);

            pertFactor = Settings.DEFAULT_JACOBIAN_PERT_FACTOR;
            for i = 1:nStates
                delta = pertFactor * (1 + abs(state(i)));

                statePlus = state;
                statePlus(i) = statePlus(i) + delta;

                stateDerivPlus = self.computeStateDeriv(statePlus);

                stateMinus = state;
                stateMinus(i) = stateMinus(i) - delta;

                stateDerivMinus = self.computeStateDeriv(stateMinus);

                A(:, i) = (stateDerivPlus - stateDerivMinus) / (2 * delta);
            end
        end


        function A = computeParamJacobian(self, state)
            nStates = self.projectile.nStates;
            nEstimatedProjectileParams = self.projectile.nEstimatedParams;
            nEstimatedPlanetParams = self.planet.nEstimatedParams;
            nEstimatedParams = nEstimatedProjectileParams + nEstimatedPlanetParams;

            A = zeros(nStates, nEstimatedParams);
            
            pertFactor = Settings.DEFAULT_JACOBIAN_PERT_FACTOR;
            for i = 1:nEstimatedProjectileParams
                paramIdx = self.projectile.estimatedParamIdxs(i);
                param = self.projectile.params(paramIdx);

                delta = pertFactor * (1 + abs(param));

                paramPlus = param + delta;
                self.projectile.params(paramIdx) = paramPlus;

                stateDerivPlus = self.computeStateDeriv(state);

                paramMinus = param - delta;
                self.projectile.params(paramIdx) = paramMinus;

                stateDerivMinus = self.computeStateDeriv(state);

                A(:, i) = (stateDerivPlus - stateDerivMinus) / (2 * delta);

                self.projectile.params(paramIdx) = param;
            end

            for i = 1:nEstimatedPlanetParams
                paramIdx = self.planet.estimatedParamIdxs(i);
                param = self.planet.params(paramIdx);

                delta = pertFactor * (1 + abs(param));

                paramPlus = param + delta;
                self.planet.params(paramIdx) = paramPlus;

                stateDerivPlus = self.computeStateDeriv(state);

                paramMinus = param - delta;
                self.planet.params(paramIdx) = paramMinus;

                stateDerivMinus = self.computeStateDeriv(state);

                A(:, nEstimatedProjectileParams + i) = (stateDerivPlus - stateDerivMinus) / (2 * delta);

                self.planet.params(paramIdx) = param;
            end
        end


        % Setters ==================================================================================

        function set.projectile(self, projectile)
            if Settings.VALIDATE_FLAG
                self.projectile = Validator.validateType(projectile, "Projectile");
            else
                self.projectile = projectile;
            end
        end

        function set.planet(self, planet)
            if Settings.VALIDATE_FLAG
                self.planet = Validator.validateType(planet, "Planet");
            else
                self.planet = planet;
            end
        end

        function set.includeParamSTM(self, includeParamSTM)
            if Settings.VALIDATE_FLAG
                self.includeParamSTM = Validator.validateType(includeParamSTM, "logical");
            else
                self.includeParamSTM = includeParamSTM;
            end
        end
    end
end