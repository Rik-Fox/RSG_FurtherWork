
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                           %
%   This is the wrapper to run Warwick HAT model on your local machine                      %
%                                                                                           %
%   Inputs:                                                                                 %
%       Cloc - a string in the format of ALPHA-3 country codes                              %
%       Ploc - an integer, provine index                                                    %
%       Zloc - an integer, health zone index                                                %
%       Aloc - an integer, health area index                                                %
%       Model - a cell array containing the model indices                                   %
%       ParaStr - a 4-digits string related to parameter setting                            %
%       MinimumData - an integer indicating the minimum requirement in Data                 %
%       RunMCMC - a boolean determining to run MCMC or not                                  %
%       RunProjection - an integer denoting the number of realizations used in Projection   %
%       RunCFS - an integer determining to run particular CounterFactual scenarios or not   %
%       RunPlot - a boolean determining to run Plot or not                                  %
%                                                                                           %    
%   Note: a Paras_xxxx.mat file stored in the same directory is required                    %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cloc = 'DRC'; % options: 'CIV', 'DRC', 'GIN', 'TCD', 'UGA'
Ploc = 1; % set as 0 if a country level run is desired
Zloc = 29; % set as 0 if a province level run is desired
Aloc = 0; % set as 0 if a heath zone level run is desired
Model = {'M4'}; % use {'M4', 'M7'} for multiple models
ParaStr = 'DRC100';
MinimumData = 10; % minimum data needed to run MCMC and Projection
RunMCMC = 0; % options: 1(run MCMC) and 0
RunProjection = 10; % options: any postive integer(number of realizations for Projection) and 0(skip) ## how many Posterior samples to use
RunSamples = 10; % options: any positive integer (number of samples per realisation, usually 10) or 0 (Projection only) ## how many stoch runs for each Posterior sample
RunReactive = 0; % options: 1(run all reactive strategies) and 0
RunCFS = 0; % options: 0(actual), 1(no vector control), 2(no mini-mobile team), 3(no Fexinidazole), 4(no RDT)
RunPlot = 1; % options: 1(plot only) and 0(simulation only)

%tic
RunLocalMachine(Cloc, Ploc, Zloc, Aloc, Model, ParaStr, MinimumData, RunMCMC, RunProjection, RunSamples, RunReactive, RunCFS, RunPlot)
%toc










