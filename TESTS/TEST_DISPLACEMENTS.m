function [uxyn1, uxyn2, pars] = TEST_DISPLACEMENTS(testnum, pars)
% This function returns a set of test displacements for asperity
% derivatives.
%
% testnum = 1-15


% a set of matching delta parameters for a usual set of parameters.
deltam = 2e-5;
deltabar = 1.692849099124275e-05;

switch testnum    
    case 1
        %Test 1 that was used for the sphere models
        uxyn1 = [0, 0, 1.9e-5]; %z must be greater than zero. 
        uxyn2 = [1e-6, 0.5e-6, 2e-5];

        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;

    case 2
        %Fully slipped test
        uxyn1 = [0, 0, 1.9e-5]; %z must be greater than zero. 
        uxyn2 = [2e-6, 5e-6, 2e-5];
                
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;

    case 3
        %Fully stuck test
        uxyn1 = [0, 0, 1.9e-5]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 2e-5];
                
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;

    case 4 
        %Fully stuck + some displacement test
        uxyn1 = [0, 0, 1.9e-5]; %z must be greater than zero. 
        uxyn2 = [0.5e-6, 0.5e-6, 2e-5];
                
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;

    case 5
        %Fully stuck test
        uxyn1 = [0, 0, 2e-5]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 2e-5-1e-6];
                
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;

    case 6
        %Fully stuck + some displacement test
        uxyn1 = [0, 0, 2e-5]; %z must be greater than zero. 
        uxyn2 = [0.5e-6, 0.5e-6, 2e-5+1e-6];
                
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;

    case 7
        %Re-establishment of contact case. - microslip
        uxyn1 = [0, 0, -2e-5]; %z must be greater than zero. 
        uxyn2 = [0.1e-6, 0.1e-6, 2e-5-1e-6];
                
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;

    case 8
        %Re-establishment of contact case. - macroslip
        uxyn1 = [0, 0, -2e-5]; %z must be greater than zero. 
        uxyn2 = [100e-6, 100e-6, 2e-5];
                
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;
        
        
    case 9
        %Fully stuck test
        uxyn1 = [0, 0, 2e-5]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 2e-5];
        disp('This case is expected to fail since z displacement is repeated leading to ill defined derivative at turning point');
                
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;

    case 10
        %Fully slipped test
        uxyn1 = [0, 0, 2e-5]; %z must be greater than zero. 
        uxyn2 = [2e-6, 5e-6, 1.9e-5];
                
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;
        
    case 11
        %Elastic Normal Regime in initial loading
        uxyn1 = [0, 0, 0]; %z must be greater than zero. 
        uxyn2 = [2e-6, 5e-6, 0.5*pars.delta_y];
                
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;
        
    case 12
        %Fully separate
        uxyn1 = [0, 0, 2e-5]; %z must be greater than zero. 
        uxyn2 = [2e-6, 5e-6, 1.0e-5];
        disp('Case is expected to give NaN''s for error, verify that the jacobian does not contain any NaN''s though')
        
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;
        
    case 13
        
        %Slight Plastic loading regime
        uxyn1 = [0, 0, pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [2e-6, 2.1e-6, 1.01*pars.delta_y];
        
        pars.mu = 0.05; %the default that most tests were written for.
        pars.sliptype = 2;
        
    case 14
        pars.mu = Inf;
        pars.sliptype = 2;
        disp('Will give Nan''s if using a model without a plasticity limit');
        
        %Fully slipped test - plasticity
        
        uxyn1 = [0, 0, 2e-5]; %z must be greater than zero. 
        uxyn2 = [2e-6, 5e-6, 1.9e-5];
        
    case 15
        pars.mu = 0.01;
        pars.sliptype = 2;
        
        %Fully slipped test - Coloumb
        uxyn1 = [0, 0, 2e-5]; %z must be greater than zero. 
        uxyn2 = [2e-6, 5e-6, 1.9e-5];
        
    case 16
        pars.mu = 0.1;
        pars.sliptype = 1;
        
        %Fully slipped test - Coloumb
        uxyn1 = [0, 0, 2e-5]; %z must be greater than zero. 
        uxyn2 = [2e-6, 5e-6, 1.9e-5];
        
    case 17
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAKE SURE TO TEST BOTH HALVES OF THE MIN IN THE CEB MODEL.
        % [0.1, 0.6, 0.99 1.2]*[ddeltam-deltabar, deltay] for second normal
        % 0, stick and slip regimes
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 0.1*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 18
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 0.6*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 19
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 0.99*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 20
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 1.2*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 21
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [1000e-6, 1000e-6, 0.1*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 22
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [1000e-6, 1000e-6, 0.6*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 23
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [1000e-6, 1000e-6, 0.99*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 24
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [1000e-6, 1000e-6, 1.2*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 25
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [.0005e-6, .001e-6, 0.1*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 26
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [0.001e-6, 0.001e-6, 0.6*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 27
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [0.001e-6, 0.001e-6, 0.99*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 28
        
        %Fully stuck test
        uxyn1 = [0, 0, 0.9*pars.delta_y]; %z must be greater than zero. 
        uxyn2 = [0.001e-6, 0.001e-6, 1.2*pars.delta_y];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
        
        
    case 29
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAKE SURE TO TEST BOTH HALVES OF THE MIN IN THE CEB MODEL.
        % [0.1, 0.6, 0.99 1.2]*[ddeltam-deltabar, deltay] for second normal
        % 0, stick and slip regimes
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 0.1*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 30
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 0.6*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 31
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 0.99*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 32
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 1.2*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 33
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [1000e-6, 1000e-6, 0.1*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 34
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [1000e-6, 1000e-6, 0.6*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 35
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [1000e-6, 1000e-6, 0.99*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 36
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [1000e-6, 1000e-6, 1.2*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 37
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [.0005e-6, .001e-6, 0.1*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 38
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [0.001e-6, 0.001e-6, 0.6*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 39
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [0.001e-6, 0.001e-6, 0.99*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 40
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [0.001e-6, 0.001e-6, 1.2*(deltam-deltabar) + deltabar];
                
        pars.mu = Inf; %the default that most tests were written for.
        pars.sliptype = 3;
        
    case 41
        
        % UTS CEB
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [0.0e-6, 0.0e-6, deltabar-2.5e-6];
                
        pars.mu = Inf; 
        pars.sliptype = 5;
        
    case 42
        
        % UTS CEB
        
        %Fully stuck test
        uxyn1 = [0, 0, deltam]; %z must be greater than zero. 
        uxyn2 = [1.0e-6, 1.0e-6, deltabar-2.5e-6];
                
        pars.mu = Inf; 
        pars.sliptype = 5;
    otherwise
        
        % Let the user know they have done all tests
        error('You have reached the maximum test number.');
        
end



end

