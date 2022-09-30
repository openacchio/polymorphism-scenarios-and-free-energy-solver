function [] = create_landscape_for_scenarios_ABC()
% [O.Penacchio, June 2022] 
% > The function generates the fitness landscapes for the paper "Outcomes of 
% multifarious selection on the evolution of visual signals",
% <add ref of paper>
% > Four scenarios (A, B1, B2 and C) are available and should be set manually
% (see "choose scenario below"); the output is a series of matrices that
% define the fitness landscape in each scenario ( 'Ztot', 'Zsex', 'Znat',
% 'Zcost'), and is saved directly to disk

% ======================================
%                   Miscellaneous
% ======================================

% ===== size/resolution of landscape =====
N = 50;
[xx, yy] = meshgrid((1:N), (1:N));

% ===== show components or not =====
showComp = 1;

% ===== choose scenario =====
% Please comment/uncomment manually to choose scenario
scenario = 'scenarioA';
scenario = 'scenarioB1';
scenario = 'scenarioB2';
scenario = 'scenarioC';

% ===== useful logistic function =====
S = @(tau, u) (1./(1 + exp(tau).^(-1*(u))));

% ======================================
%    MAIN:    get corresponding landscape 
% ======================================

% ===== Zcost: Free energy for phenotype cost =====
% N.B.1: The cost is the same across scenarios and therefore defined here
% once for all.
% N.B.2: For the cost, we need to give a prohibitive value for phenotypes that correspond
% to a high values for both colour and luminance contrast as these phenotypes 
% ***are not feasible***! (See Crothers and Cumming 2013 for details; all measurements
% follow and "antidiagona"l of the form Col_contrast + Lum_contrast = constant, which seems to
% be the way the effect of the colorations is maximised)
beta = 3.0;
kappa = 2.25;
Zcost = kappa*( 1.4* ((xx/N) + (yy/N) )/2 ).^(beta);
Zcost(Zcost > 3) = 3;
Zcost = Zcost  +  0.75*( fliplr( (xx/N).^9 ) + flipud( (yy/N).^9 ));

% ===== get other components (Znat and Zsex) =====
switch scenario
    case 'scenarioA'
        % ===== Zsex: Free energy for sexual selection =====
        lambdaSex = 1;
        deltaSex = 0.55; 
        thetaSex = 90 + 37.5; 
        Z0 = flipud( S(0.25, cosd(thetaSex)*(xx- (N - deltaSex*N)) + sind(thetaSex)*(yy - (N - deltaSex*N))) );
        Zsex = lambdaSex*Z0;
        
        % ===== Znat: Free energy for natural selection =====
        lambdaNat = 1;
        Znat = lambdaNat*Z0';
        
        % ===== plot components and output =====
        if showComp == 1
            figure('Position', [50, 50, 1200, 850]),
            subplot(2,3,1),surf(xx, yy, Zsex), view([0, 90]), title('Zsex')
            subplot(2,3,2),surf(xx, yy, Znat), view([0, 90]), title('Znat')
            subplot(2,3,3),surf(xx, yy, Zsex + Znat), view([0, 90]), title('Zsex + Znat')
            subplot(2,3,4),surf(xx, yy, Zcost), view([0, 90]), title('Zcost')
            Ztot = Zsex + Znat + Zcost;
            subplot(2,3,6),surf(xx, yy, Ztot), view([0, 90]), title('full landscape (Ztot)')
        end
        
    case 'scenarioB1'
        % ===== Zsex: Free energy for sexual selection =====
        lambdaSex = 0.5;
        deltaSex = 0.55;
        thetaSex = 90;
        Z0 = flipud( S(0.25, cosd(thetaSex)*(xx- (N - deltaSex*N)) + sind(thetaSex)*(yy - (N - deltaSex*N))) );
        Zsex = lambdaSex*Z0;
        
        % ===== Znat: Free energy for natural selection =====
        lambdaNat = 1.5;
        Znat = lambdaNat*Z0';
                       
        % ===== plot components and output =====
        if showComp == 1
            figure('Position', [50, 50, 1200, 850]),
            subplot(2,3,1),surf(xx, yy, Zsex), view([0, 90]), title('Zsex')
            subplot(2,3,2),surf(xx, yy, Znat), view([0, 90]), title('Znat')
            subplot(2,3,3),surf(xx, yy, Zsex + Znat), view([0, 90]), title('Zsex + Znat')
            subplot(2,3,4),surf(xx, yy, Zcost), view([0, 90]), title('Zcost')
            Ztot = Zsex + Znat + Zcost;
            subplot(2,3,6),surf(xx, yy, Ztot), view([0, 90]), title('full landscape (Ztot)')
        end
        
    case 'scenarioB2'
        % ===== Zsex: Free energy for sexual selection =====
        lambdaSex = 1.5;
        deltaSex = 0.55; 
        thetaSex = 90; 
        Z0 = flipud( S(0.25, cosd(thetaSex)*(xx- (N - deltaSex*N)) + sind(thetaSex)*(yy - (N - deltaSex*N))) );
        Zsex = lambdaSex*Z0;
        
        % ===== Znat: Free energy for natural selection =====
        lambdaNat = 0.5;
        Znat = lambdaNat*Z0';
                        
        % ===== plot components and output =====
        if showComp == 1
            figure('Position', [50, 50, 1200, 850]),
            subplot(2,3,1),surf(xx, yy, Zsex), view([0, 90]), title('Zsex')
            subplot(2,3,2),surf(xx, yy, Znat), view([0, 90]), title('Znat')
            subplot(2,3,3),surf(xx, yy, Zsex + Znat), view([0, 90]), title('Zsex + Znat')
            subplot(2,3,4),surf(xx, yy, Zcost), view([0, 90]), title('Zcost')
            Ztot = Zsex + Znat + Zcost;
            subplot(2,3,6),surf(xx, yy, Ztot), view([0, 90]), title('full landscape (Ztot)')
        end
        
    case 'scenarioC'
        % ===== Zsex: Free energy for sexual selection =====
        lambdaSex = 1;
        deltaSex = 0.55; % ref: 0.55
        thetaSex = 90; % ref 55
        Z0 = flipud( S(0.25, cosd(thetaSex)*(xx- (N - deltaSex*N)) + sind(thetaSex)*(yy - (N - deltaSex*N))) );
        Zsex = lambdaSex*Z0;
        
        % ===== Znat: Free energy for natural selection =====
        lambdaNat = 1;
        Znat = lambdaNat*Z0';
        
        % ===== plot components and output =====
        if showComp == 1
            figure('Position', [50, 50, 1200, 850]),
            subplot(2,3,1),surf(xx, yy, Zsex), view([0, 90]), title('Zsex')
            subplot(2,3,2),surf(xx, yy, Znat), view([0, 90]), title('Znat')
            subplot(2,3,3),surf(xx, yy, Zsex + Znat), view([0, 90]), title('Zsex + Znat')
            subplot(2,3,4),surf(xx, yy, Zcost), view([0, 90]), title('Zcost')
            Ztot = Zsex + Znat + Zcost;
            subplot(2,3,6),surf(xx, yy, Ztot), view([0, 90]), title('full landscape (Ztot)')
        end
        
end

% ===========================================
%       save output for Python's PDE treatment
% ===========================================
save([pwd '\total_free_energy_landscape_for_Fokker_Plank_equation_' scenario '_Siz' num2str(N) '.mat'], 'Ztot', 'Zsex', 'Znat', 'Zcost', '-v7')

end