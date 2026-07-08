classdef Mars < Planet
    properties (Constant)
        DEFAULT_G = 3.728;
        
        DEFAULT_RHO = 0.01332;
        DEFAULT_RHO0 = 0.01332;
        DEFAULT_RHO_TABLE_X = [0; 25000];
        DEFAULT_RHO_TABLE_Y = [0.01332; 0.01332];

        % DEFAULT_R = 188.9;
        % DEFAULT_T = 223;
        % DEFAULT_GAMMA = 1.335 
        DEFAULT_A = 237.1;  % TODO: Compute using R and T as params instead

        DEFAULT_H = 12547;       % TODO: Scaling (rescale to km?)
        
        DEFAULT_VWINDX = 0;
        DEFAULT_VWINDX_TABLE_X = [0; 25000];
        DEFAULT_VWINDX_TABLE_Y = [0; 0];
        DEFAULT_VWINDY = 0;
        DEFAULT_VWINDY_TABLE_X = [0; 25000];
        DEFAULT_VWINDY_TABLE_Y = [0; 0];
    end


    methods
        % Constructor ==============================================================================

        function self = Mars(varargin)
            self = self@Planet(varargin{:});
        end
    end
end