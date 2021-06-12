function out = def(action,varargin)

codedir = fileparts(mfilename('fullpath'));
simdir = fullfile(codedir,'simulations');
if ~exist(simdir,'dir'), mkdir(simdir); end

switch action
    case 'pipedir'
        out = simdir;        
    case 'addpath'
        addpath(fullfile(codedir,'tools'));        
        
    case 'unp'
        out = 'Unpredictability';
    case 'vol'
        out = 'Volatility';
    case 'lr'
        out = 'Learning rate';
    case 'col1'
        out = [1 .2 .2];
    case 'col'
        out = [0    0.4470    0.7410; 0.8500    0.3250    0.0980];
%         out = [.2 0 0; .8 .4 .1; 1 .2 .2] ;
%         out = [.8 .4 .1; 1 .2 .2] ;

    case 'col_br'
        col = ones(1,3)*.1;
%         col(2,:) = [255 85 155]/255;
        col(2,:) = [255 85 85]/255;      
        col(2,:) = [215 25 25]/255;
            out = col;

    case 'alf'
        out = .6;
    case 'fs'
        out = 12;
    case 'fsy'
        out = 14;
    case 'fsA'
        out = 18;
    case 'fn'
        out = 'Arial';
    case 'xsA'        
        out = -.2;
    case 'ysA'
        out = 1.05;
    case 'abc'
        out = 'abcdefghijklmnopqrstuv';
    
    otherwise
        error('%s does not exist',action);
end
end

