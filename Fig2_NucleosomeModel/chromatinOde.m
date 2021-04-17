% This file is automatically generated by update_model.m from the chromatin_model spreadsheet.
% Model URL: F:/enhancer_dynamics/model_v2/chromatin_model.xlsx

%% Section 1: Declaration/initialization (code from ode_init.m)
function delta = chromatinOde(t,x,time,v, tf)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% delta = chromatinOde(t,x,ode_options,v)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Full chromatin ODE model - stiff system, designed to be solved using ode15 (see this function's help file to set
% options).  
% In phase 1 (v.PHASE == 1),simulation is run w/o stimulus, until convergence to initialize steady-state,
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set the perisistant variable that discrete-delay rxns refer to 

% persistent DELAY;
% if isempty(t)
%     % The sim function calls the ode function once without arguments to reset persistent variables
%     sz = 5000; % Size of total delay memory - if DELAYS.idx exceeds this, function will throw error
%     DELAY.t = zeros(sz,1);
%     DELAY.E0 = zeros(sz,1); % increase if getting errors
%     DELAY.idx     = 1; % index (starts at 1)
%     return;
% end


% Slice parameters, get previous concentrations
[v.PARAMS, v.SPECIES] = chromatinInitialize();
p = v.PARAMS;
delta = zeros(size(x));

tf = interp1(time, tf, t, 'linear');



%% Section 2: Unpack species
E_0 = x(1);
E_1 = x(2);
E_2 = x(3);
E_3 = x(4);
E_4 = x(5);
E_5 = x(6);
E_6 = x(7);
E_7 = x(8);
E_8 = x(9);
E_9 = x(10);
E_10 = x(11);
E_11 = x(12);
E_12 = x(13);
E_13 = x(14);
E_14 = x(15);


%% Section 4: Set reaction rates
rxn_1 = p(1,1) * (tf.^p(1,2))/( (tf.^p(1,2)) + (p(1,3).^p(1,2)) ) * E_0;
rxn_2 = p(2,1) * E_1;
rxn_3 = p(3,1) * (tf.^p(3,2))/( (tf.^p(3,2)) + (p(3,3).^p(3,2)) ) * E_1;
rxn_4 = p(4,1) * E_2;
rxn_5 = p(5,1) * (tf.^p(5,2))/( (tf.^p(5,2)) + (p(5,3).^p(5,2)) ) * E_2;
rxn_6 = p(6,1) * E_3;
rxn_7 = p(7,1) * (tf.^p(7,2))/( (tf.^p(7,2)) + (p(7,3).^p(7,2)) ) * E_3;
rxn_8 = p(8,1) * E_4;
rxn_9 = p(9,1) * (tf.^p(9,2))/( (tf.^p(9,2)) + (p(9,3).^p(9,2)) ) * E_4;
rxn_10 = p(10,1) * E_5;
rxn_11 = p(11,1) * (tf.^p(11,2))/( (tf.^p(11,2)) + (p(11,3).^p(11,2)) ) * E_5;
rxn_12 = p(12,1) * E_6;
rxn_13 = p(13,1) * (tf.^p(13,2))/( (tf.^p(13,2)) + (p(13,3).^p(13,2)) ) * E_6;
rxn_14 = p(14,1) * E_7;
rxn_15 = p(15,1) * (tf.^p(15,2))/( (tf.^p(15,2)) + (p(15,3).^p(15,2)) ) * E_7;
rxn_16 = p(16,1) * E_8;
rxn_17 = p(17,1) * (tf.^p(17,2))/( (tf.^p(17,2)) + (p(17,3).^p(17,2)) ) * E_8;
rxn_18 = p(18,1) * E_9;
rxn_19 = p(19,1) * (tf.^p(19,2))/( (tf.^p(19,2)) + (p(19,3).^p(19,2)) ) * E_9;
rxn_20 = p(20,1) * E_10;
rxn_21 = p(21,1) * (tf.^p(21,2))/( (tf.^p(21,2)) + (p(21,3).^p(21,2)) ) * E_10;
rxn_22 = p(22,1) * E_11;
rxn_23 = p(23,1) * (tf.^p(23,2))/( (tf.^p(23,2)) + (p(23,3).^p(23,2)) ) * E_11;
rxn_24 = p(24,1) * E_12;
rxn_25 = p(25,1) * (tf.^p(25,2))/( (tf.^p(25,2)) + (p(25,3).^p(25,2)) ) * E_12;
rxn_26 = p(26,1) * E_13;
rxn_27 = p(27,1) * (tf.^p(27,2))/( (tf.^p(27,2)) + (p(27,3).^p(27,2)) ) * E_13;
rxn_28 = p(28,1) * E_14;


%% Section 5: Set species' deltas from reactions
delta(1) = - rxn_1 + rxn_2;
delta(2) = - rxn_2 - rxn_3 + rxn_1 + rxn_4;
delta(3) = - rxn_4 - rxn_5 + rxn_3 + rxn_6;
delta(4) = - rxn_6 - rxn_7 + rxn_5 + rxn_8;
delta(5) = - rxn_8 - rxn_9 + rxn_7 + rxn_10;
delta(6) = - rxn_10 - rxn_11 + rxn_9 + rxn_12;
delta(7) = - rxn_12 - rxn_13 + rxn_11 + rxn_14;
delta(8) = - rxn_14 - rxn_15 + rxn_13 + rxn_16;
delta(9) = - rxn_16 - rxn_17 + rxn_15 + rxn_18;
delta(10) = - rxn_18 - rxn_19 + rxn_17 + rxn_20;
delta(11) = - rxn_20 - rxn_21 + rxn_19 + rxn_22;
delta(12) = - rxn_22 - rxn_23 + rxn_21 + rxn_24;
delta(13) = - rxn_24 - rxn_25 + rxn_23 + rxn_26;
delta(14) = - rxn_26 - rxn_27 + rxn_25 + rxn_28;
delta(15) = - rxn_28 + rxn_27;
