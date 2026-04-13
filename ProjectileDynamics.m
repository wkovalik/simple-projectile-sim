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
        function [Fx, Fy] = projectileGravModel(self)
            m = self.projectile.params.m.value;

            g = self.earth.params.g.value;
            
            Fx = 0;
            Fy = -m * g;
        end

        function [Fx, Fy] = projectileAeroModel(self, state)          
            vx = state(3);
            vy = state(4);
            
            S = self.projectile.params.S.value;
            CD = self.projectile.params.CD.value;

            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;
            
            vInfx = vx - vWindx;
            VInf = (vInfx ^ 2 + vy ^ 2) ^ 0.5;
            
            Fx = -(rho * S * CD / 2) * VInf * vInfx;
            Fy = -(rho * S * CD / 2) * VInf * vy;
        end
        
        % Derivative methods =======================================================================
        function stateDeriv = computeStateDeriv(self, state)
            vx = state(3);
            vy = state(4);
            
            m = self.projectile.params.m.value;
            
            [FGravx, FGravy] = self.projectileGravModel();
            [FAerox, FAeroy] = self.projectileAeroModel(state);
            
            Fx = FGravx + FAerox;
            Fy = FGravy + FAeroy;
            
            xDeriv = vx;
            yDeriv = vy;
            vxDeriv = Fx / m;
            vyDeriv = Fy / m;
            
            stateDeriv = [xDeriv; yDeriv; vxDeriv; vyDeriv];
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
            vx = state(3);
            vy = state(4);

            m = self.projectile.params.m.value;
            S = self.projectile.params.S.value;
            CD = self.projectile.params.CD.value;

            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;

            vInfx = vx - vWindx;
            VInf = (vInfx ^ 2 + vy ^ 2) ^ 0.5;

            A = zeros(self.N_PROJECTILE_STATES);

            A(1, 3) = 1;
            A(3, 3) = -(rho * S * CD / (2 * m)) * (vInfx ^ 2 / VInf + VInf);
            A(4, 3) = -(rho * S * CD / (2 * m)) * (vInfx * vy / VInf);

            A(2, 4) = 1;
            A(3, 4) = -(rho * S * CD / (2 * m)) * (vInfx * vy / VInf);
            A(4, 4) = -(rho * S * CD / (2 * m)) * (vy ^ 2 / VInf + VInf);
        end

        function A = computeParamJacobian(self, state)
            A = zeros(self.N_PROJECTILE_STATES, self.N_ESTIMATED_PARAMS);

            for i = 1:self.N_ESTIMATED_PARAMS
                A(:, i) = self.estimatedParamJacobianMap{i}(self, state);  % See Note 1
            end
        end

        function A = computeDragJacobian(self, state)
            vx = state(3);
            vy = state(4);
            
            m = self.projectile.params.m.value;
            S = self.projectile.params.S.value;

            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;
            
            vInfx = vx - vWindx;
            VInf = (vInfx ^ 2 + vy ^ 2) ^ 0.5;
            
            A = zeros(self.N_PROJECTILE_STATES, 1);
            
            A(3, 1) = -(rho * S / (2 * m)) * VInf * vInfx;
            A(4, 1) = -(rho * S / (2 * m)) * VInf * vy;
        end
            
        function A = computeWindJacobian(self, state)
            vx = state(3);
            vy = state(4);
            
            m = self.projectile.params.m.value;
            S = self.projectile.params.S.value;
            CD = self.projectile.params.CD.value;

            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;
            
            vInfx = vx - vWindx;
            VInf = (vInfx ^ 2 + vy ^ 2) ^ 0.5;
            
            A = zeros(self.N_PROJECTILE_STATES, 1);
            
            A(3, 1) = (rho * S * CD / (2 * m)) * (vInfx ^ 2 / VInf + VInf);
            A(4, 1) = (rho * S * CD / (2 * m)) * (vInfx * vy / VInf);
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