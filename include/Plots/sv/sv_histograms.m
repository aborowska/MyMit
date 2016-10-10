 figure(600)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
	set(gcf,'defaulttextinterpreter','latex');
    suptitle('adapt')
    subplot(3,1,1)
    hist(theta_adapt(:,1),20)
    title('$$c$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(3,1,2)
    hist(theta_adapt(:,2),20)
    title('$$\phi$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(3,1,3)
    hist(theta_adapt(:,3),20)
    title('$$\sigma^{2}_{\eta}$$')
	set(gca,'TickLabelInterpreter','latex')
    
    
 figure(700)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
	set(gcf,'defaulttextinterpreter','latex');
    suptitle('new')
    subplot(3,1,1)
    hist(theta_new(:,1),20)
    title('$$c$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(3,1,2)
    hist(theta_new(:,2),20)
    title('$$\phi$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(3,1,3)
    hist(theta_new(:,3),20)
    title('$$\sigma^{2}_{\eta}$$')
	set(gca,'TickLabelInterpreter','latex')    
    
figure(200)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
	set(gcf,'defaulttextinterpreter','latex');
    suptitle('1')
    subplot(3,1,1)
    hist(theta1(:,1),20)
    title('$$c$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(3,1,2)
    hist(theta1(:,2),20)
    title('$$\phi$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(3,1,3)
    hist(theta1(:,3),20)
    title('$$\sigma^{2}_{\eta}$$')
	set(gca,'TickLabelInterpreter','latex')
    
    
figure(300)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
	set(gcf,'defaulttextinterpreter','latex');
    suptitle('hl')
    subplot(3,1,1)
    hist(theta_hl(:,1),20)
    title('$$c$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(3,1,2)
    hist(theta_hl(:,2),20)
    title('$$\phi$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(3,1,3)
    hist(theta_hl(:,3),20)
    title('$$\sigma^{2}_{\eta}$$')
	set(gca,'TickLabelInterpreter','latex')