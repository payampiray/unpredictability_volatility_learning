function data = data_amygdala_costa(all_data)
if nargin<1, all_data = 0; end

% Table S1 of the paper

% 60/40
% Control
% 0.68 (Acq), 0.67 (Rev)
% Amygdala
% 0.57 (Acq), 0.54 (Rev)
% 
% 100/0
% Control
% 0.91 (Acq), 0.94 (Rev)
% Amygdala
% 0.77 (Acq), 0.66 (Rev)

specs=[     {'Control-Stochastic'}    {'Control-Deterministic'}    {'Lesioned-Stochastic'}    {'Lesioned-Deterministic'}
            {'Control'           }    {'Control'              }    {'Lesioned'           }    {'Lesioned'              }
            {'Stochastic'        }    {'Deterministic'        }    {'Stochastic'         }    {'Deterministic'         }
            {'60/40'             }    {'100/0'                }    {'60/40'              }    {'100/0'                 }
            {[                 1]}    {[                    2]}    {[                  3]}    {[                     4]}];


x_acq =    [0.68    0.91    0.57    0.77];
x_rev =    [0.67    0.94    0.54    0.66];

e_acq = 0.01*ones(1,4);
e_rev = 0.01*ones(1,4);

data = struct('x_acq',x_acq,'x_rev',x_rev,'e_acq',e_acq,'e_rev',e_rev,'specs',{specs});

if all_data
    data = full_data;
end

end

function data = full_data

x_acq = [.68 .77 .86 .91 .57 .62 .66 .77];
x_rev = [.67 .81 .86 .94 .54 .57 .61 .66];

specs = [ {'control-60','control-70','control-80','control-100', 'Lesioned-60','Lesioned-70','Lesioned-80','Lesioned-100'} ];

e_acq = 0.01*ones(1,8);
e_rev = 0.01*ones(1,8);

data = struct('x_acq',x_acq,'x_rev',x_rev,'e_acq',e_acq,'e_rev',e_rev,'specs',{specs});
end