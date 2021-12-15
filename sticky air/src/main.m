% planetesimal: main model routine

% print run header
fprintf(1,'\n\n************************************************************\n');
fprintf(1,    '*****  planetesimal  |  %s  |  %s  *****\n'         ,RUN.ID,datetime);
fprintf(1,    '************************************************************\n\n');


% %% initialise model run
% initialise;


%% physical time stepping loop

while NUM.time <= NUM.tend && NUM.step <= NUM.maxstep
    
    % print time step header
    fprintf(1,'\n*****  step = %d;  dt = %1.4e;  time = %1.4e yr \n\n',NUM.step,NUM.dt/NUM.yr,NUM.time/NUM.yr);
    
    % % =============================================================
    % % Solve fluid mechanics equations
    % % =============================================================
    
    if ~mod(NUM.step,round(2*RUN.nup/NUM.CFL))  % perform every 'nup' gridsteps of transport
        if RUN.selfgrav; solve_gravity; end     % update gravity
        solve_fluidmech;                        % solve  fluid mechanics
        up2date;
        
    end
    
    %     figure(4)
    %     subplot(1,3,1)
    %     imagesc(NUM.xP,NUM.zP,DEF.txx); colormap(subplot(1,3,1), cm2); colorbar;
    %     axis ij image; title('x normal  stress')
    %
    %     subplot(1,3,2)
    %     imagesc(NUM.xP,NUM.zP,DEF.tzz); colormap(subplot(1,3,2), cm2); colorbar;
    %     axis ij image; title('z normal stress')
    %
    %     subplot(1,3,3)
    %     imagesc(NUM.xP,NUM.zP,DEF.txz); colormap(subplot(1,3,3), cm2); colorbar;
    %     axis ij image; title('xz shear stress')
    
    %update level set free surface
    Div_PHI0 = Div_PHI;
    Div_PHI = advection(SOL.U,SOL.W,NUM.PHI,NUM.dx,NUM.dz,'third upwind');
%     NUM.PHI(2:end-1,2:end-1) = NUM.PHI(2:end-1,2:end-1) + Div_PHI.*NUM.dt;
    NUM.PHI(2:end-1,2:end-1) = NUM.PHI(2:end-1,2:end-1) - (NUM.theta .*Div_PHI   ...
        + (1-NUM.theta).*    Div_PHI0) .* NUM.dt;
    % apply boundary conditions
    % apply top boundary conditions  
    NUM.PHI(1,:) = NUM.PHI (2,:);
    
    % apply bottom boundary conditions
    NUM.PHI(end,:) = NUM.PHI(end-1,:);
    
    % apply side boundary conditions
    NUM.PHI(:,[1 end]) = NUM.PHI(:,[2 end-1]);
    
    DPHI = full(Div_PHI);
    %     figure(3)
    %     imagesc(NUM.xP,NUM.zP,NUM.PHI); axis ij image; colormap(cm2), colorbar
    % % =============================================================
    % % Solve thermo-chemical equations
    % % =============================================================
    
    solve_thermochem;               % solve every time step
    
    
    % % =============================================================
    % % Update material & deformation properties
    % % =============================================================
    
    up2date;                                % update materials & deformation
    
    
    if ~mod(NUM.step,round(2*RUN.nop/NUM.CFL))
        output;                     % output every 'nop' gridsteps transport
%         figure(3)
%         imagesc(NUM.xP,NUM.zP,SOL.P); axis ij image; colormap(cm2), colorbar
%         title('pressure')
%         drawnow
    end
    
    % increment time
    NUM.step = NUM.step + 1;
    NUM.time = NUM.time + NUM.dt;
end