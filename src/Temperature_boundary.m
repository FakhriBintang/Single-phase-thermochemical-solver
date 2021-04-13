% update temperature boundary

% interpolate temperature onto staggered grid
T_vx(:,2:nx1) = (T_mid(:,2:nx1)+T_mid(:,1:nx))./2;
T_vz(2:nz1,:) = (T_mid(2:nz1,:)+T_mid(1:nz,:))./2;

% Apply thermal boundary condition to interpolated nodes
% upper boundary, constant temperature
T_vx(1,:)    = 2*T_top   - T_vx(2,:);
T_vz(1,:)   = 2*T_top   - T_vz(2,:);

% lower boundary, constant temperature
T_vx(end,:)   = 2*(T_bot+300)   - T_vx(nz1,:);
T_vz(end,:) = 2*(T_bot+300)   - T_vz(nz,:);
%       T_vx(end,:) = T_bot;
%       T_vz(end,:) = T_bot;

% left boundary, insulating boundary
T_vx(:,1)       = T_vx(:,2);
T_vz(1:nz1,1)   = T_vz(1:nz1,2);
% right boundary, insulating boundary
T_vx(:,nx1)     = T_vx(:,nx);
T_vz(1:nz1,Nx)  = T_vz(1:nz1,nx1);
switch BotBoundary
    case 'hell portal'
        % plume beneath lithosphere
        midx = fix(L/2);
        indexP = find(abs(xp(1,:)) <= midx+ (L*0.25) & abs(xp(1,:)) >= midx - (L*0.25));
        % indexP = 1:Nx; indexvx = 1:Nx; indexvz = 1:Nx;
        T_mid(end-1:end,indexP) = T_bot+300;
        indexvx = find(abs(xvx(1,:)) <= midx+ (L*0.25) & abs(xvx(1,:)) >= midx - (L*0.25));
        T_vx(end-1:end,indexvx) = T_bot+300;
        indexvz = find(abs(xvz(1,:)) <= midx+ (L*0.25) & abs(xvz(1,:)) >= midx - (L*0.25));
        T_vz(end-1:end,indexvz) = T_bot+300;
        
    case 'uniform'
end
[Rho_vx,Rho_vz,Rho_mid] = Thermal_expansion(Rho0,...
    T_mid,T_vx,T_vz,Alpha_vx,Alpha_vz,Alpha_mid,T_top);