classdef ProjectileDynamics < handle
    properties
        projectile
        earth

        estimatedParams
        estimatedParamSourceMap
        estimatedParamAnalyticJacobianMap

        jacobianMethod
    end

    properties (Dependent)
        N_PROJECTILE_STATES
        N_ESTIMATED_PARAMS
    end

    properties (Constant)
        PARAM_SOURCE_MAP = dictionary( ...
                                 "CD",     "projectile", ...
                                 "vWindx", "earth",      ...
                                 "vWindy", "earth"       ...
                             );
        PARAM_ANALYTIC_JACOBIAN_MAP = dictionary( ...
                                          "CD",     @computeAnalyticDragJacobian,  ...
                                          "vWindx", @computeAnalyticWindxJacobian, ...
                                          "vWindy", @computeAnalyticWindyJacobian  ...
                                      );
    end

    methods
        % Constructor ==============================================================================
        function self = ProjectileDynamics(projectile, earth)
            % Set handles for projectile and planet
            self.projectile = projectile;
            self.earth = earth;
            
            % Initialize list of parameters to be estimated (used to compute parameter Jacobian)
            self.estimatedParams = [];
            
            % Set initial Jacobian computation method
            self.jacobianMethod = Constants.DEFAULT_JACOBIAN_METHOD;
        end


        % Force model methods ======================================================================
        function [Fx, Fy, Fz] = projectileGravModel(self)
            % Get projectile parameters
            m = self.projectile.params.m.value;
            
            % Get planet parameters
            g = self.earth.params.g.value;
            
            % Compute force vector components
            Fx = 0;
            Fy = 0;
            Fz = m * g;
        end

        function [Fx, Fy, Fz] = projectileAeroModel(self, state)
            % Get projectile state
            vx = state(4);
            vy = state(5);
            vz = state(6);
            
            % Get projectile parameters
            S = self.projectile.params.S.value;
            CD = self.projectile.params.CD.value;
            
            % Get planet parameters
            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;
            vWindy = self.earth.params.vWindy.value;
            
            % Compute atmosphere-relative velocities
            vInfx = vx - vWindx;
            vInfy = vy - vWindy;
            vInfz = vz;
            VInf = (vInfx ^ 2 + vInfy ^ 2 + vInfz ^ 2) ^ 0.5;
            
            % Compute force vector components
            Fx = -(rho * S * CD / 2) * VInf * vInfx;
            Fy = -(rho * S * CD / 2) * VInf * vInfy;
            Fz = -(rho * S * CD / 2) * VInf * vInfz;
        end
        

        % Derivative methods =======================================================================
        function stateDeriv = computeStateDeriv(self, state)
            % Get projectile state
            vx = state(4);
            vy = state(5);
            vz = state(6);
            
            % Get projectile parameters
            m = self.projectile.params.m.value;
            
            % Compute forces
            [FGravx, FGravy, FGravz] = self.projectileGravModel();
            [FAerox, FAeroy, FAeroz] = self.projectileAeroModel(state);
            
            Fx = FGravx + FAerox;
            Fy = FGravy + FAeroy;
            Fz = FGravz + FAeroz;
            
            % Compute state derivatives
            xDeriv = vx;
            yDeriv = vy;
            zDeriv = vz;
            vxDeriv = Fx / m;
            vyDeriv = Fy / m;
            vzDeriv = Fz / m;
            
            stateDeriv = [xDeriv; yDeriv; zDeriv; vxDeriv; vyDeriv; vzDeriv];
        end
        
        function [stateSTMDeriv, paramSTMDeriv] = computeSTMDeriv(self, state, stateSTM, paramSTM)
            % Note: stateSTM = d(projectile states)/d(initial projectile states)
            %       paramSTM = d(projectile states)/d(estimated parameters)
            %
            %       stateA = d(projectile state derivatives)/d(projectile states)
            %       paramA = d(projectile state derivatives)/d(estimated parameters)
            
            % Compute state STM derivative ---------------------------------------------------------

            % Unvectorize STM
            stateSTM = reshape(stateSTM, [self.N_PROJECTILE_STATES, self.N_PROJECTILE_STATES]);
            
            % Compute dynamics Jacobian
            switch self.jacobianMethod
                case "numeric"
                    stateA = self.computeNumericStateJacobian(state);
                case "analytic"
                    stateA = self.computeAnalyticStateJacobian(state);
            end
            
            % Compute STM derivative
            stateSTMDeriv = stateA * stateSTM;

            % Revectorize STM
            stateSTMDeriv = stateSTMDeriv(:);
            
            % Compute parameter STM derivative -----------------------------------------------------

            if ~isempty(self.estimatedParams)
                % Unvectorize STM
                paramSTM = reshape(paramSTM, [self.N_PROJECTILE_STATES, self.N_ESTIMATED_PARAMS]);
                
                % Compute dynamics Jacobian
                switch self.jacobianMethod
                    case "numeric"
                        paramA = self.computeNumericParamJacobian(state);
                    case "analytic"
                        paramA = self.computeAnalyticParamJacobian(state);
                end
                    
                % Compute STM derivative
                paramSTMDeriv = stateA * paramSTM + paramA;

                % Revectorize STM
                paramSTMDeriv = paramSTMDeriv(:);

            else
                paramSTMDeriv = [];

            end
        end

        function augStateDeriv = computeStateAndSTMDeriv(self, augState)
            % Deconcatenate augmented state
            state = augState(1:self.N_PROJECTILE_STATES);
            stateSTM = augState(self.N_PROJECTILE_STATES + (1:self.N_PROJECTILE_STATES ^ 2));
            if ~isempty(self.estimatedParams)
                paramSTM = augState((self.N_PROJECTILE_STATES + self.N_PROJECTILE_STATES ^ 2) + (1:(self.N_PROJECTILE_STATES * self.N_ESTIMATED_PARAMS)));
            else
                paramSTM = [];
            end
            
            % Compute projectile state derivative
            stateDeriv = self.computeStateDeriv(state);

            % Compute vectorized STM derivatives
            [stateSTMDeriv, paramSTMDeriv] = self.computeSTMDeriv(state, stateSTM, paramSTM);
            
            % Reconcatenate into augmented state
            augStateDeriv = [stateDeriv; stateSTMDeriv; paramSTMDeriv];
        end
        

        % Jacobian methods =========================================================================
        function A = computeNumericStateJacobian(self, state)
            A = zeros(self.N_PROJECTILE_STATES);
            
            % Numerically compute partial derivatives w.r.t. each projectile state
            for i = 1:self.N_PROJECTILE_STATES
                % Compute state perturbation
                delta = Constants.JACOBIAN_PERTURBATION_FACTOR * (1 + abs(state(i)));
                
                % Compute forward-perturbed state derivative
                statePlus = state;
                statePlus(i) = statePlus(i) + delta;

                stateDerivPlus = self.computeStateDeriv(statePlus);
                
                % Compute backward-perturbed state derivative
                stateMinus = state;
                stateMinus(i) = stateMinus(i) - delta;

                stateDerivMinus = self.computeStateDeriv(stateMinus);
                
                % Approximate partial derivative using central difference method
                A(:, i) = (stateDerivPlus - stateDerivMinus) / (2 * delta);
            end
        end

        function A = computeAnalyticStateJacobian(self, state)
            % Get projectile states
            vx = state(4);
            vy = state(5);
            vz = state(6);

            % Get projectile parameters
            m = self.projectile.params.m.value;
            S = self.projectile.params.S.value;
            CD = self.projectile.params.CD.value;

            % Get planet parameters
            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;
            vWindy = self.earth.params.vWindy.value;

            % Compute atmosphere-relative velocities
            vInfx = vx - vWindx;
            vInfy = vy - vWindy;
            vInfz = vz;
            VInf = (vInfx ^ 2 + vInfy ^ 2 + vInfz ^ 2) ^ 0.5;

            % Build Jacobian -----------------------------------------------------------------------
            A = zeros(self.N_PROJECTILE_STATES);

            % Partial derivatives w.r.t. vx
            A(1, 4) = 1;
            A(4, 4) = -(rho * S * CD / (2 * m)) * (vInfx ^ 2 / VInf + VInf);
            A(5, 4) = -(rho * S * CD / (2 * m)) * (vInfx * vInfy / VInf);
            A(6, 4) = -(rho * S * CD / (2 * m)) * (vInfx * vInfz / VInf);

            % Partial derivatives w.r.t. vy
            A(2, 5) = 1;
            A(4, 5) = -(rho * S * CD / (2 * m)) * (vInfx * vInfy / VInf);
            A(5, 5) = -(rho * S * CD / (2 * m)) * (vInfy ^ 2 / VInf + VInf);
            A(6, 5) = -(rho * S * CD / (2 * m)) * (vInfy * vInfz / VInf);

            % Partial derivatives w.r.t. vz
            A(3, 6) = 1;
            A(4, 6) = -(rho * S * CD / (2 * m)) * (vInfx * vInfz / VInf);
            A(5, 6) = -(rho * S * CD / (2 * m)) * (vInfy * vInfz / VInf);
            A(6, 6) = -(rho * S * CD / (2 * m)) * (vInfz ^ 2 / VInf + VInf);
        end

        function A = computeNumericParamJacobian(self, state)
            A = zeros(self.N_PROJECTILE_STATES, self.N_ESTIMATED_PARAMS);

            % Numerically compute partial derivatives w.r.t. each estimated parameter
            for i = 1:self.N_ESTIMATED_PARAMS
                paramName = self.estimatedParams{i};
                paramSource = self.estimatedParamSourceMap{i};

                param = self.(paramSource).params.(paramName).value;

                % Compute parameter perturbation
                delta = Constants.JACOBIAN_PERTURBATION_FACTOR * (1 + param);

                % Compute forward-perturbed state derivative
                paramPlus = param + delta;
                self.(paramSource).params.(paramName).value = paramPlus;

                stateDerivPlus = self.computeStateDeriv(state);

                % Compute backward-perturbed state derivative
                paramMinus = param - delta;
                self.(paramSource).params.(paramName).value = paramMinus;

                stateDerivMinus = self.computeStateDeriv(state);

                % Approximate partial derivative using central difference method
                A(:, i) = (stateDerivPlus - stateDerivMinus) / (2 * delta);

                % Reset parameter
                self.(paramSource).params.(paramName).value = param;
            end
        end

        function A = computeAnalyticParamJacobian(self, state)
            A = zeros(self.N_PROJECTILE_STATES, self.N_ESTIMATED_PARAMS);

            % Build Jacobian using individual Jacobians corresponding to each estimated parameter
            for i = 1:self.N_ESTIMATED_PARAMS
                A(:, i) = self.estimatedParamAnalyticJacobianMap{i}(self, state);  % See Note 1
            end
        end

        function A = computeAnalyticDragJacobian(self, state)
            % Get projectile state
            vx = state(4);
            vy = state(5);
            vz = state(6);

            % Get projectile parameters
            m = self.projectile.params.m.value;
            S = self.projectile.params.S.value;

            % Get planet parameters
            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;
            vWindy = self.earth.params.vWindy.value;

            % Compute atmosphere-relative velocities
            vInfx = vx - vWindx;
            vInfy = vy - vWindy;
            vInfz = vz;
            VInf = (vInfx ^ 2 + vInfy ^ 2 + vInfz ^ 2) ^ 0.5;

            % Build Jacobian -----------------------------------------------------------------------
            A = zeros(self.N_PROJECTILE_STATES, 1);

            % Partial derivatives w.r.t. CD
            A(4, 1) = -(rho * S / (2 * m)) * VInf * vInfx;
            A(5, 1) = -(rho * S / (2 * m)) * VInf * vInfy;
            A(6, 1) = -(rho * S / (2 * m)) * VInf * vInfz;
        end
            
        function A = computeAnalyticWindxJacobian(self, state)
            % Get projectile state
            vx = state(4);
            vy = state(5);
            vz = state(6);

            % Get projectile parameters
            m = self.projectile.params.m.value;
            S = self.projectile.params.S.value;
            CD = self.projectile.params.CD.value;

            % Get planet parameters
            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;
            vWindy = self.earth.params.vWindy.value;

            % Compute atmosphere-relative velocities
            vInfx = vx - vWindx;
            vInfy = vy - vWindy;
            vInfz = vz;
            VInf = (vInfx ^ 2 + vInfy ^ 2 + vInfz ^ 2) ^ 0.5;

            % Build Jacobian -----------------------------------------------------------------------
            A = zeros(self.N_PROJECTILE_STATES, 1);

            % Partial derivatives w.r.t. vWindx
            A(4, 1) = (rho * S * CD / (2 * m)) * (vInfx ^ 2 / VInf + VInf);
            A(5, 1) = (rho * S * CD / (2 * m)) * (vInfx * vInfy / VInf);
            A(6, 1) = (rho * S * CD / (2 * m)) * (vInfx * vInfz / VInf);
        end

        function A = computeAnalyticWindyJacobian(self, state)
            % Get projectile state
            vx = state(4);
            vy = state(5);
            vz = state(6);

            % Get projectile parameters
            m = self.projectile.params.m.value;
            S = self.projectile.params.S.value;
            CD = self.projectile.params.CD.value;

            % Get planet parameters
            rho = self.earth.params.rho.value;
            vWindx = self.earth.params.vWindx.value;
            vWindy = self.earth.params.vWindy.value;

            % Compute atmosphere-relative velocities
            vInfx = vx - vWindx;
            vInfy = vy - vWindy;
            vInfz = vz;
            VInf = (vInfx ^ 2 + vInfy ^ 2 + vInfz ^ 2) ^ 0.5;

            % Build Jacobian -----------------------------------------------------------------------
            A = zeros(self.N_PROJECTILE_STATES, 1);

            % Partial derivatives w.r.t. vWindy
            A(4, 1) = (rho * S * CD / (2 * m)) * (vInfx * vInfy / VInf);
            A(5, 1) = (rho * S * CD / (2 * m)) * (vInfy ^ 2 / VInf + VInf);
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

            % Rebuild mapping lists for all estimated parameters
            self.buildEstimatedParamMaps();
        end

        function buildEstimatedParamMaps(self)
            % Builds list of object source location and list of Jacobian compute function handles
            % for each corresponding estimated parameter. See Note 1
            % Example: self.estimatedParams = [ "CD", "vWindx" ]
            %          self.estimatedParamSourceMap = { "projectile", "earth" }
            %          self.estimatedParamJacobianMap = { @self.computeDragJacobian, @self.computeWindxJacobian }

            self.estimatedParamSourceMap = {};
            self.estimatedParamAnalyticJacobianMap = {};

            if ~isempty(self.estimatedParams)
                for i = 1:self.N_ESTIMATED_PARAMS
                    estimatedParam = self.estimatedParams(i);
                    
                    self.estimatedParamSourceMap{end + 1} = self.PARAM_SOURCE_MAP(estimatedParam);
                    self.estimatedParamAnalyticJacobianMap{end + 1} = self.PARAM_ANALYTIC_JACOBIAN_MAP(estimatedParam);
                end
            end
        end

        function set.jacobianMethod(self, jacobianMethod)
            self.jacobianMethod = Validator.validateInEnum(jacobianMethod, ["analytic", "numeric"]);
        end
    end
end

% TODO Tomorrow
% Move maps to here. Add param source map. Then compute param jacobian like state jacobian
% Do same for sensormodel. Will have to add projectile and earth handles to sensormodels
% Try implementing TODO about not passing state vectors, but rather projectile (which contains state) and earth models?
%     see if performance takes a huge hit
% TODO Seems doing param assigment is very slow. Maybe concatenate all params into vector (like
% state) and pass that around (so all funcs are now fn(state, params) or fn(projState, projParams, planetParams) instead of fn(state))


% Note 1
%
% The original implementation of the computeAnalyticParamJacobian() function was
%
% | for i = 1:self.N_ESTIMATED_PARAMS
% |     estimatedParam = self.estimatedParams(i);
% | 
% |     A(:, i) = feval(self.PARAM_ANALYTIC_JACOBIAN_MAP(estimatedParam), state);
% | end
%
% However, directly looking up the Jacobian function handles in self.PARAM_ANALYTIC_JACOBIAN_MAP
% ended up being quite slow (around twice as slow as the current implementation, in fact). I dunno
% why. My best guess is the local, individual self.computeAnalyticParamJacobian function handles
% were being rebuilt on every single lookup of self.PARAM_ANALYTIC_JACOBIAN_MAP(estimatedParam).
%
% So instead, a local self.estimatedParamAnalyticJacobianMap property is built, which is a cell
% array of the subset of Jacobian function handles in the full self.PARAM_ANALYTIC_JACOBIAN_MAP for
% the parameters in self.estimatedParams only. This cell array lists the handles in the correct
% order according to the order in self.estimatedParams and is rebuilt whenever self.estimatedParams
% is updated.
%
% For example:
% 
% | self.estimatedParams = ["CD", "vWindx"]
% | self.estimatedParamAnalyticJacobianMap = { @computeAnalyticDragJacobian, @computeAnalyticWindxJacobian }
%
% This way, the function handles in self.estimatedParamAnalyticJacobianMap are already fully defined
% when calling self.computeAnalyticParamJacobian(), and since self.estimatedParamAnalyticJacobianMap
% is a cell array, it is simply iterated over. Similar reasoning applies for building
% self.estimatedParamSourceMap.