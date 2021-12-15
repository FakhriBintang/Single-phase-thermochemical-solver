% solve_flag
%% P-cells 
tic
NUM.flag.P = zeros(NUM.nzP,NUM.nxP);
% deactivate cells above the free surface
if RUN.selfgrav
    NUM.flag.P(NUM.RP>NUM.D*0.9/2) = 3;
else
    NUM.flag.P(NUM.RP<NUM.D/10) = 3;
end
% activate Free surface cells (need to figure out a way to vextorize)

for j = 2:1:NUM.nxP-1
for i = 2:1:NUM.nzP-1
    
    if    (NUM.flag.P(i,j) == 3 && NUM.flag.P(i+1,j)  == 0 ...
        || NUM.flag.P(i,j) == 3 && NUM.flag.P(i-1,j)  == 0 ...
        || NUM.flag.P(i,j) == 3 && NUM.flag.P(i,j-1)  == 0 ...
        || NUM.flag.P(i,j) == 3 && NUM.flag.P(i,j+1)  == 0)
        
        NUM.flag.P(i,j) = 6;
    end
end
end


%% Vx-nodes; 0 = active 
NUM.flag.U = ones(NUM.nzU,NUM.nxU);
NUM.flag.U = NUM.flag.P(:,1:end-1).*NUM.flag.P(:,2:end);
NUM.flag.U(NUM.flag.U>0) = 3;

%% Vz-nodes; 0 = active
NUM.flag.W = zeros(NUM.nzW,NUM.nxW);
NUM.flag.W = NUM.flag.P(1:end-1,:).*NUM.flag.P(2:end,:);
NUM.flag.W(NUM.flag.W>0) = 3;

%% corner nodes
NUM.flag.C = zeros(NUM.nzC,NUM.nxC);
NUM.flag.C = (NUM.flag.P(1:end-1,1:end-1) ...
    +  NUM.flag.P(2:end  ,1:end-1) ...
    +  NUM.flag.P(1:end-1,2:end  ) ...
    +  NUM.flag.P(2:end  ,2:end  ));
NUM.flag.C(NUM.flag.C>0) = 3;
toc
figure(10)
subplot(2,2,1)
imagesc(NUM.xP,NUM.zP,NUM.flag.P)
colorbar
% hold on
% plot(NUM.xP,NUM.Fs,'--k', 'linewidth', 3)
title('P-flagging')
axis tight

subplot(2,2,2)
imagesc(NUM.xU,NUM.zU,NUM.flag.U)
colorbar
% hold on
% plot(NUM.xP,NUM.Fs,'--k', 'linewidth', 3)
title('U-flagging')
axis tight

subplot(2,2,3)
imagesc(NUM.xW,NUM.zW,NUM.flag.W)
colorbar
% hold on
% plot(NUM.xP,NUM.Fs,'--k', 'linewidth', 3)
title('W-flagging')
axis tight

subplot(2,2,4)
imagesc(NUM.xC,NUM.zC,NUM.flag.C)
colorbar
% hold on
% plot(NUM.xP,NUM.Fs,'--k', 'linewidth', 3)
title('C-flagging')
axis tight
