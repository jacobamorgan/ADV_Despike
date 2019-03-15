% close all; clear variables; clc
% fldr = 'Z:\JAM\Flume\F01_R03\velocity\';%20161115\X06m\Y40cm\';
% fls = dir([fldr '*.mat']);
% for i = 1:1%size(fls,1)
%     load([fldr fls(i).name])
%     V = Data.Profiles_VelX;
% end
function Vclean = adv_despike(V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program name: adv_despike
% Author: Jacob A. Morgan
%         Colorado State University
%         Department of Civil & Environmental Engineering
% Date: 10 October 2017
% Description: This program consists of a function implementing the
%     despiking of acoustic Doppler velocimeter data as presented in Goring
%     and Nikora (2002) and modified by Wahl (2003). The function is set up
%     to handle the input of velocity profile data of dimensions M x N,
%     where M is the length of the time-series and N is the number of
%     vertical nodes. This program does not consider the correlation or
%     signal-to-noise ratio of the input data. Currently the function
%     replaces detected spikes by peforming a cubic interpolation of 12
%     valid velocity data points on either side of the spike.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vworking = V;
iis = 1:5;%5;
for ii = iis
    disp(['iteration ' num2str(ii) ' of ' num2str(iis(end))])
    if ii > 1
        Vworking = Vclean;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 0: Initialize velocity data, filter out long-term trend(s)
    meanV = mean(Vworking);
    windowSize = 100*5; % 100 Hz, windowSize is filter window
    v2 = filter(ones(1,round(windowSize))/round(windowSize),1,...
        Vworking-repmat(meanV,size(Vworking,1),1));
    newV = Vworking-v2;
    medianV = median(newV);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 1: Calculate surrogates for the first and second derivatives.
    U(:,:,1) = 100*(newV-repmat(medianV,size(newV,1),1)); % use median (not mean), per Wahl (2003)
    for i = 2:3
        for j = 1:size(U,2)
            U(:,j,i) = gradient(U(:,j,i-1));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2: a) Calculate the standard deviations of all three variables.
    sU = 1.483*median(abs(U-repmat(median(U),size(U,1),1))); % use MAD, per Wahl (2003)
    %         b) Calculate the expected maximum using Universal criterion
    lamU = (2*log(size(U,1)))^0.5; % Universal threshold
    maxU = squeeze(lamU.*sU);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 3: Calculate the rotation angle of the principal axis of d2u
    %         versus u using the cross correlation.
    alph = atan(sum(U(:,:,1).*U(:,:,3))./sum(U(:,:,1).^2)); % theta in GN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 4: Calculate the ellipse that has maxima and minima from above.
    aa = [1 2 1]; bb = [2 3 3];
    sig = repmat(0:2*pi/149:2*pi,size(maxU,1),1);
    for i = 1:3
        a = aa(i); b = bb(i);
        r = (repmat(maxU(:,a).^2.*maxU(:,b).^2,1,size(sig,2))./...
            (repmat(maxU(:,a).^2,1,size(sig,2)).*(sin(sig)).^2+...
            repmat(maxU(:,b).^2,1,size(sig,2)).*(cos(sig)).^2)).^0.5;
        [xi,yi] = pol2cart(sig,r);
        if i == 3 % rotate ellipse in u-d2u axis by alpha
            x(:,:,i) = xi.*repmat(cos(alph'),1,size(xi,2))-...
                yi.*repmat(sin(alph'),1,size(yi,2));
            y(:,:,i) = xi.*repmat(sin(alph'),1,size(xi,2))+...
                yi.*repmat(cos(alph'),1,size(yi,2));
        else
            x(:,:,i) = xi; y(:,:,i) = yi;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 5: Identify points that lie outside ellipse and replace them.
    %         a) Convert into spherical coordinates, per Wahl (2003), in
    %            order to identify location of points in reference to
    %            ellipsoid. This makes everything a lot faster, and maybe
    %            more accurate too (?).
    mU = mean(U(:,:,1));
    U(:,:,1) = U(:,:,1)-repmat(mU,size(U,1),1); % center velocity pts around zero
    rhU = sqrt(U(:,:,1).^2+U(:,:,2).^2+U(:,:,3).^2); % rho of velocity pts
    thU = atan(U(:,:,2)./U(:,:,1)); % theta of velocity pts
    thU(U(:,:,1)<0) = pi+thU(U(:,:,1)<0); % quadrant is important for theta
    phU = atan(sqrt(U(:,:,1).^2+U(:,:,2).^2)./U(:,:,3)); % phi velocity pts
    rh2in = (sin(phU).*cos(thU).*repmat(cos(alph),size(phU,1),1)+...
        cos(phU).*repmat(sin(alph),size(phU,1),1)).^2./...
        repmat((maxU(:,1).^2)',size(phU,1),1)+...
        (sin(phU).*cos(thU).*repmat(sin(alph),size(phU,1),1)-...
        cos(phU).*repmat(cos(alph),size(phU,1),1)).^2./...
        repmat((maxU(:,2).^2)',size(phU,1),1)+...
        (sin(phU).*sin(thU)).^2./...
        repmat((maxU(:,3).^2)',size(phU,1),1); % inv sq rh of ellipsoid
    rhe = 1./sqrt(rh2in); % rho of ellipsoid at velocity angles
    ine = rhU<=rhe; % if true, velocity point falls inside ellipsoid
    U(:,:,1) = U(:,:,1)+repmat(mU,size(U,1),1); % restore velocity pts
    %         b) Replace velocity "spike" points with "corrected" points
    Vclean = U(:,:,1); % create new velocity variable to manipulate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%        IS THERE A BETTER/FASTER WAY TO DO ALL THIS??        %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    replacemethod = 1;
    if replacemethod == 1
        for z = 1:size(U,2)
            indzout = find(ine(:,z)==0);
            indzin = find(ine(:,z)==1);
            for j = 1:length(indzout)
                k = find(indzin<=indzout(j),1,'last');
                if indzout(j)<indzin(1); k = 1; end
                if k > 11
                    tbefore = indzin(k-11:k);
                    vbefore = Vclean(indzin(k-11:k),z);
                else
                    tbefore = indzin(1:k);
                    vbefore = Vclean(indzin(1:k),z);
                end
                if k < length(indzin)-12
                    tafter = indzin(k+1:k+12);
                    vafter = Vclean(indzin(k+1:k+12),z);
                else
                    if k == length(indzin)
                        tafter = []; vafter = [];
                    else
                        tafter = indzin(k+1:end);
                        vafter = Vclean(indzin(k+1:end),z);
                    end
                end
                vspike = Vclean(indzout(j),z);
                vreplace = interp1([tbefore;tafter],[vbefore;vafter],...
                    indzout(j),'pchip'); % perform cubic interpolation
                if abs(vreplace) > abs(vspike)
                    vreplace = vspike;
                end
                Vclean(indzout(j),z) = vreplace;
            end
        end
    end
    Vclean = Vclean./100+repmat(medianV,size(v2,1),1)+v2;
end
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % z = 20;
% z = 35;
% inez = ine(:,z);
% fx = 500;
% fy = 500;
% set(0,'Units','points'); ss = get(0,'screensize');
% figure('Units','points',...
%     'Position',[(ss(3)-fx)/2 (ss(4)-fy)/2 fx fy],...
%     'PaperUnits','points','PaperSize',[fx fy],...
%     'PaperPosition',[0 0 fx fy],...
%     'DefaultAxesFontSize',12);
% s = [1 4 3];
% for i = 1:3
%     a = aa(i); b = bb(i);
%     sp(i) = subplot(2,2,s(i)); hold on; box on
%     plot(squeeze(U(inez,z,a)),squeeze(U(inez,z,b)),'.k','MarkerSize',4)
%     plot(squeeze(U(~inez,z,a)),squeeze(U(~inez,z,b)),'.r','MarkerSize',4)
%     plot(squeeze(x(z,:,i)),squeeze(y(z,:,i)),'r','LineWidth',1.5)
% end
% set(sp(3),'Position',get(sp(3),'Position')+[0.02 0.02 0.04 0.06])
% sp3p = get(sp(3),'Position');
% set(sp(1),'Position',[sp3p(1) sum(sp3p([2 4])) sp3p(3) sp3p(4)])
% set(sp(2),'Position',[sum(sp3p([1 3])) sp3p(2) sp3p(3) sp3p(4)])
% set(sp,'TickLength',[0.02 0.02])
% for i = 1:3
%     if i == 1
%         set(sp(i),'XTickLabel','')
%         ytl = get(sp(i),'YTickLabel'); ytl{1} = '';
%         set(sp(i),'YTickLabel',ytl);
%         ylabel(sp(i),'\Delta{\itu} (cm/s)')
%     elseif i == 2
%         set(sp(i),'YTickLabel','')
%         xtl = get(sp(i),'XTickLabel'); xtl{1} = '';
%         set(sp(i),'XTickLabel',xtl);
%         xlabel(sp(i),'\Delta{\itu} (cm/s)')
%     elseif i == 3
%         xlabel('{\itu} (cm/s)')
%         ylabel(sp(i),'\Delta^2{\itu} (cm/s)')
%     end
% end
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t = 0:0.01:359.99;
% z = 18;
% inez = ine(:,z);
% fx = 1200;
% fy = 300;
% set(0,'Units','points'); ss = get(0,'screensize');
% figure('Units','points',...
%     'Position',[(ss(3)-fx)/2 (ss(4)-fy)/2 fx fy],...
%     'PaperUnits','points','PaperSize',[fx fy],...
%     'PaperPosition',[0 0 fx fy],...
%     'DefaultAxesFontSize',12);
% hold on; box on
% plot(t,V(:,z),'r')
% plot(t,Vclean(:,z),'k')
% % plot(t,smooth(V(:,z),6000),'g')
% % plot(t,Vclean(:,z),'.c')
% plot(t(~inez),V(~inez,z),'.c','MarkerSize',10)
% plot(t(~inez),Vclean(~inez,z),'.g','MarkerSize',10)
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t = 0:0.01:359.99;
% z = 20;
% inez = ine(:,z);
% fx = 900;
% fy = 500;
% set(0,'Units','points'); ss = get(0,'screensize');
% figure('Units','points',...
%     'Position',[(ss(3)-fx)/2 (ss(4)-fy)/2 fx fy],...
%     'PaperUnits','points','PaperSize',[fx fy],...
%     'PaperPosition',[0 0 fx fy],...
%     'DefaultAxesFontSize',12);
% subplot(2,1,1)
% hold on; box on
% plot(t,V(:,z),'r')
% ylm = get(gca,'YLim');
% subplot(2,1,2)
% hold on; box on
% plot(t,Vclean(:,z),'k','MarkerSize',10)
% set(gca,'YLim',ylm)