% Planetesimal: main model routine

% print run header
fprintf(1,'\n\n************************************************************\n');
fprintf(1,    '*****  planetesimal  |  %s  |  %s  *****\n'         ,RUN.ID,datetime);
fprintf(1,    '************************************************************\n\n');


%% initialise model run
initialise;


%% physical time stepping loop

while NUM.time < NUM.tend && NUM.step < NUM.maxstep

    % print time step header
    fprintf(1,'\n*****  step = %d;  time = %1.4e yr \n\n',NUM.step,NUM.time/NUM.yr);
    
    % % =============================================================
    % % Solve fluid mechanics equations
    % % =============================================================
    
    if ~mod(NUM.step,round(2*RUN.nup/NUM.CFL))
        solve_fluidmech;            % solve  every 'nup' gridsteps of transport
        up2date;                    % update every 'nup' gridsteps of transport
    end
    
    % % =============================================================
    % % Solve thermo-chemical equations
    % % =============================================================
    
    solve_thermochem;               % solve every time step
    
    if ~mod(NUM.step,round(2*RUN.nop/NUM.CFL)) 
        output;                     % output every 'nop' gridsteps transport
    end
    
    % increment time
    NUM.step = NUM.step + 1;
    NUM.time = NUM.time + NUM.dt;
end