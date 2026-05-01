classdef ProjectileDynamics < handle
    % TODO: Add tempStorage to derivative functions. That way not recomputing a ton of variables on each force model or Jacobian call

    properties
        projectile
        planet

        includeParamSTM = false;
    end
    

    methods
        % Constructor ==============================================================================

        function self = ProjectileDynamics(projectile, planet)
            if nargin ~= 2
                error("Not enough input arguments. Requires projectile and planet.")
            end

            self.projectile = projectile;
            self.planet = planet;
        end

        
        % Force model methods ======================================================================

        function F = computeGravityForce(self)
            m = self.projectile.params(self.projectile.mIdx);

            g = self.planet.computeGravity();

            F = zeros(3, 1);
            F(3) = m * g;
        end


        function F = computeAeroForce(self, state)
            z = state(3);
            vx = state(4);
            vy = state(5);
            vz = state(6);

            S = self.projectile.params(self.projectile.SIdx);
            
            h = -z;
            [rho, a] = self.planet.computeAtmosphere(h);
            [vWindx, vWindy] = self.planet.computeWind(h);
            
            vAtmx = vx - vWindx;
            vAtmy = vy - vWindy;
            vAtmz = vz;
            VAtm = (vAtmx ^ 2 + vAtmy ^ 2 + vAtmz ^ 2) ^ 0.5;

            mach = VAtm / a;

            CD = self.projectile.computeAeroCoeffs(mach);
            
            k = -(rho * S * CD / 2);

            F = zeros(3, 1);
            F(1) = k * VAtm * vAtmx;
            F(2) = k * VAtm * vAtmy;
            F(3) = k * VAtm * vAtmz;
        end

        
        % Derivative methods =======================================================================

        function stateDeriv = computeStateDeriv(self, state)
            nStates = self.projectile.nStates;

            vx = state(4);
            vy = state(5);
            vz = state(6);
            
            % TODO: Compute acceleration function?
            m = self.projectile.params(self.projectile.mIdx);

            FGrav = self.computeGravityForce();
            FAero = self.computeAeroForce(state);

            F = FGrav + FAero;

            stateDeriv = zeros(nStates, 1);
            stateDeriv(1) = vx;
            stateDeriv(2) = vy;
            stateDeriv(3) = vz;
            stateDeriv(4) = F(1) / m;
            stateDeriv(5) = F(2) / m;
            stateDeriv(6) = F(3) / m;
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

        function A = computeStateJacobian(self, state)
            nStates = self.projectile.nStates;

            A = zeros(nStates);
            
            pertFactor = 1E-04;
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
            
            pertFactor = 1E-04;
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