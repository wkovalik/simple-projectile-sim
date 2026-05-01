classdef Earth < Planet
    properties (Constant)
        DEFAULT_G = 9.81;
        
        DEFAULT_RHO = 1.225;
        DEFAULT_RHO0 = 1.225;
        DEFAULT_H = 8500;       % TODO: Scaling (rescale to km?)
        % DEFAULT_R = 287.05287;
        % DEFAULT_T = 288.15;
        DEFAULT_A = 344.1;  % TODO: Compute using R and T as params instead

        DEFAULT_VWINDX = 0;
        DEFAULT_VWINDX_TABLE_X = [0; 10000];
        DEFAULT_VWINDX_TABLE_Y = [0; 0];
        DEFAULT_VWINDY = 0;
        DEFAULT_VWINDY_TABLE_X = [0; 10000];
        DEFAULT_VWINDY_TABLE_Y = [0; 0];
    end


    methods
        % Constructor ==============================================================================

        function self = Earth(varargin)
            self = self@Planet(varargin{:});
        end
    end
end