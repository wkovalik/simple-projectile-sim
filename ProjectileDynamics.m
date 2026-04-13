classdef ProjectileDynamics < handle
    properties
        projectile
        earth

        estimatedParams
        estimatedParamJacobianMap
    end

    properties (Dependent)
        N_PROJECTILE_STATES
        N_ESTIMATED_PARAMS
    end

    methods
        % Constructor ==============================================================================
        function self = ProjectileDynamics(projectile, earth)
            self.projectile = projectile;
            self.earth = earth;

            self.estimatedParams = [];
        end

        % Force model methods ======================================================================
        function [Fx, Fy, Fz] = projectileGravModel(self)
            m = self.projectile.params.m.value;

            g = self.earth.params.g.value;
            
            Fx = 0;
            Fy = 0;
            Fz = m * g;
        end

        function [Fx, Fy, Fz] = projectileAeroModel(self, state)          
            vx = state(4);
            vy = state(5);
            vz = state(6);
            
            S = self.projectile.params.S.value;
            CD = self.projectile.params.CD.value;

            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;
            
            vInfx = vx - vWindx;
            vInfy = vy;
            vInfz = vz;
            VInf = (vInfx ^ 2 + vInfy ^ 2 + vInfz ^ 2) ^ 0.5;
            
            Fx = -(rho * S * CD / 2) * VInf * vInfx;
            Fy = -(rho * S * CD / 2) * VInf * vInfy;
            Fz = -(rho * S * CD / 2) * VInf * vInfz;
        end
        
        % Derivative methods =======================================================================
        function stateDeriv = computeStateDeriv(self, state)
            vx = state(4);
            vy = state(5);
            vz = state(6);
            
            m = self.projectile.params.m.value;
            
            [FGravx, FGravy, FGravz] = self.projectileGravModel();
            [FAerox, FAeroy, FAeroz] = self.projectileAeroModel(state);
            
            Fx = FGravx + FAerox;
            Fy = FGravy + FAeroy;
            Fz = FGravz + FAeroz;
            
            xDeriv = vx;
            yDeriv = vy;
            zDeriv = vz;
            vxDeriv = Fx / m;
            vyDeriv = Fy / m;
            vzDeriv = Fz / m;
            
            stateDeriv = [xDeriv; yDeriv; zDeriv; vxDeriv; vyDeriv; vzDeriv];
        end
        
        function [stateSTMDeriv, paramSTMDeriv] = computeSTMDeriv(self, state, stateSTM, paramSTM)
            stateSTM = reshape(stateSTM, [self.N_PROJECTILE_STATES, self.N_PROJECTILE_STATES]);

            stateA = self.computeStateJacobian(state);

            stateSTMDeriv = stateA * stateSTM;
            stateSTMDeriv = stateSTMDeriv(:);

            if ~isempty(self.estimatedParams)
                paramSTM = reshape(paramSTM, [self.N_PROJECTILE_STATES, self.N_ESTIMATED_PARAMS]);
                
                paramA = self.computeParamJacobian(state);
    
                paramSTMDeriv = stateA * paramSTM + paramA;
                paramSTMDeriv = paramSTMDeriv(:);

            else
                paramSTMDeriv = [];

            end
        end

        function augStateDeriv = computeStateAndSTMDeriv(self, augState)
            state = augState(1:self.N_PROJECTILE_STATES);
            stateSTM = augState(self.N_PROJECTILE_STATES + (1:self.N_PROJECTILE_STATES ^ 2));
            if ~isempty(self.estimatedParams)
                paramSTM = augState((self.N_PROJECTILE_STATES + self.N_PROJECTILE_STATES ^ 2) + (1:(self.N_PROJECTILE_STATES * self.N_ESTIMATED_PARAMS)));
            else
                paramSTM = [];
            end

            stateDeriv = self.computeStateDeriv(state);
            [stateSTMDeriv, paramSTMDeriv] = self.computeSTMDeriv(state, stateSTM, paramSTM);

            augStateDeriv = [stateDeriv; stateSTMDeriv; paramSTMDeriv];
        end
        
        % Jacobian methods =========================================================================
        function A = computeStateJacobian(self, state)
            vx = state(4);
            vy = state(5);
            vz = state(6);

            m = self.projectile.params.m.value;
            S = self.projectile.params.S.value;
            CD = self.projectile.params.CD.value;

            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;

            vInfx = vx - vWindx;
            vInfy = vy;
            vInfz = vz;
            VInf = (vInfx ^ 2 + vInfy ^ 2 + vInfz ^ 2) ^ 0.5;

            A = zeros(self.N_PROJECTILE_STATES);

            A(1, 4) = 1;
            A(4, 4) = -(rho * S * CD / (2 * m)) * (vInfx ^ 2 / VInf + VInf);
            A(5, 4) = -(rho * S * CD / (2 * m)) * (vInfx * vInfy / VInf);
            A(6, 4) = -(rho * S * CD / (2 * m)) * (vInfx * vInfz / VInf);

            A(2, 5) = 1;
            A(4, 5) = -(rho * S * CD / (2 * m)) * (vInfx * vInfy / VInf);
            A(5, 5) = -(rho * S * CD / (2 * m)) * (vInfy ^ 2 / VInf + VInf);
            A(6, 5) = -(rho * S * CD / (2 * m)) * (vInfy * vInfz / VInf);

            A(3, 6) = 1;
            A(4, 6) = -(rho * S * CD / (2 * m)) * (vInfx * vInfz / VInf);
            A(5, 6) = -(rho * S * CD / (2 * m)) * (vInfy * vInfz / VInf);
            A(6, 6) = -(rho * S * CD / (2 * m)) * (vInfz ^ 2 / VInf + VInf);
        end

        function A = computeParamJacobian(self, state)
            A = zeros(self.N_PROJECTILE_STATES, self.N_ESTIMATED_PARAMS);

            for i = 1:self.N_ESTIMATED_PARAMS
                A(:, i) = self.estimatedParamJacobianMap{i}(self, state);  % See Note 1
            end
        end

        function A = computeDragJacobian(self, state)
            vx = state(4);
            vy = state(5);
            vz = state(6);
            
            m = self.projectile.params.m.value;
            S = self.projectile.params.S.value;

            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;
            
            vInfx = vx - vWindx;
            vInfy = vy;
            vInfz = vz;
            VInf = (vInfx ^ 2 + vInfy ^ 2 + vInfz ^ 2) ^ 0.5;
            
            A = zeros(self.N_PROJECTILE_STATES, 1);
            
            A(4, 1) = -(rho * S / (2 * m)) * VInf * vInfx;
            A(5, 1) = -(rho * S / (2 * m)) * VInf * vInfy;
            A(6, 1) = -(rho * S / (2 * m)) * VInf * vInfz;
        end
            
        function A = computeWindJacobian(self, state)
            vx = state(4);
            vy = state(5);
            vz = state(6);
            
            m = self.projectile.params.m.value;
            S = self.projectile.params.S.value;
            CD = self.projectile.params.CD.value;

            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;
            
            vInfx = vx - vWindx;
            vInfy = vy;
            vInfz = vz;
            VInf = (vInfx ^ 2 + vInfy ^ 2 + vInfz ^ 2) ^ 0.5;
            
            A = zeros(self.N_PROJECTILE_STATES, 1);
            
            A(4, 1) = (rho * S * CD / (2 * m)) * (vInfx ^ 2 / VInf + VInf);
            A(5, 1) = (rho * S * CD / (2 * m)) * (vInfx * vInfy / VInf);
            A(6, 1) = (rho * S * CD / (2 * m)) * (vInfx * vInfz / VInf);
        end

        % Getters ==================================================================================
        function N_PROJECTILE_STATES = get.N_PROJECTILE_STATES(self)
            N_PROJECTILE_STATES = self.projectile.state.N_STATES;
        end

        function N_ESTIMATED_PARAMS = get.N_ESTIMATED_PARAMS(self)
            N_ESTIMATED_PARAMS = length(self.estimatedParams);
        end
        
        % Setters ==================================================================================
        function set.projectile(self, projectile)
            self.projectile = Validator.validateType(projectile, "Projectile");
        end

        function set.earth(self, earth)
            self.earth = Validator.validateType(earth, "Earth");
        end

        function set.estimatedParams(self, estimatedParams)
            if ~isempty(estimatedParams)
                self.estimatedParams = Validator.validateType(estimatedParams, "string");
            end

            self.updateEstimatedParamJacobianMap();
        end

        function updateEstimatedParamJacobianMap(self)
            self.estimatedParamJacobianMap = {};

            if ~isempty(self.estimatedParams)
                for i = 1:self.N_ESTIMATED_PARAMS
                    estimatedParam = self.estimatedParams(i);
    
                    self.estimatedParamJacobianMap{end + 1} = Constants.PARAM_JACOBIAN_MAP(estimatedParam);
                end
            end
        end
    end
end

% Note 1
%
% The original implementation of the computeParamJacobian() function was
%
% | for i = 1:self.N_ESTIMATED_PARAMS
% |     estimatedParam = self.estimatedParams(i);
% | 
% |     A(:, i) = feval(Constants.PARAM_JACOBIAN_MAP(estimatedParam), self, state);
% | end
%
% However, directly looking up the Jacobian function handles in Constants.PARAM_JACOBIAN_MAP ended
% up being quite slow (around twice as slow as the current implementation, in fact). I dunno why.
% My best guess is the local, individual self.computeEstimatedParamJacobian function handles were
% being rebuilt on every single lookup of Constants.PARAM_JACOBIAN_MAP(estimatedParam).
%
% So instead, a local self.estimatedParamJacobianMap property is built which is a cell array of the
% subset of Jacobian function handles in the full Constants.PARAM_JACOBIAN_MAP for the parameters in
% self.estimatedParams only. This cell array lists the handles in the correct order according to the
% order in self.estimatedParams and is rebuilt whenever self.estimatedParams is updated.
%
% For example:
% 
% | self.estimatedParams = ["CD", "vWindx"]
% | self.estimatedParamJacobianMap = { @computeDragJacobian, @computeWindJacobian }
%
% This way, the function handles in self.estimatedParamJacobianMap are already fully defined when
% calling self.computeParamJacobian(), and since self.estimatedParamJacobianMap is a cell array, it
% is simply iterated over.