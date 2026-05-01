classdef Settings
    properties (Constant)
        DEFAULT_TIME_TOL = 1E-09;  % See Note 1

        DEFAULT_HISTORY_LEN = 1E+05;

        DEFAULT_INTERP_METHOD = "linear";
        DEFAULT_EXTRAP_METHOD = "linear";  % See Note 2
        
        DEFAULT_STATE_JACOBIAN_METHOD = "analytic";
        DEFAULT_PARAM_JACOBIAN_METHOD = "numeric";
        DEFAULT_JACOBIAN_PERT_FACTOR = 1E-04;

        DEFAULT_MAX_ITERS = 10;
        DEFAULT_CONVERGENCE_TOL = 1E-04;

        VALIDATE_FLAG = false;  % Keep on while developing. Turn off for faster performance.
    end
end


% Note 1:
%
% Subtracting this tolerance prevents any timing issues due to imprecise equality or inequality
% checks using floating point numbers.
% 
% For example, the main propagation loop checks if time < finalTime, and it should ideally terminate
% when time >= finalTime. Suppose finalTime = 7.0 and time = 6.9999999999997. This really should be
% the final integration step and the propagation loop should end here. But since time < finalTime,
% the propagator erroneously steps one more time. Using time < (finalTime - TIME_TOLERANCE) fixes
% this issue, where TIME_TOLERANCE << Integrator.stepPeriod.

% Note 2:
%
% This prevents NaN when interpolating at initial time and/or final time due to floating point
% imprecision thinking these are erroneously outside range of interpolation times.
% 
% TODO: Add warning if this ever occurs
