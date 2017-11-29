%{
11/28/2017
FFT for growth curve
relative standard deviation for different groups
%}


%% Extract data
% strains raw growth curves
B=xlsread('strainB.xlsx');
C=xlsread('strainC.xlsx');
D=xlsread('strainD.xlsx');
E=xlsread('strainE.xlsx');
F=xlsread('strainF.xlsx');
G=xlsread('strainG.xlsx');
H=xlsread('strainH.xlsx');
t=B(:,13);
%intergrate growth curve
Gro_Curv=[];
for i=1 : 12
    Gro_Curv(:,i)=B(:,i);
    Gro_Curv(:,(i+12))=C(:,i);
    Gro_Curv(:,(i+12*2))=D(:,i);
    Gro_Curv(:,(i+12*3))=E(:,i);
    Gro_Curv(:,(i+12*4))=F(:,i);
    Gro_Curv(:,(i+12*5))=G(:,i);
    Gro_Curv(:,(i+12*6))=H(:,i);
end

% 7 strains raw growth rate
for i=1:(7*12)
    RGR(:,i)=gradient(Gro_Curv(:,i))./gradient(t);
end
t_diff=t(2:length(t));

%% plot raw grwoth rate and growth curve
% figure(1);
% for i=1:7
%     subplot(3,7,i);
%     for j=(1+(i-1)*12):(i*12)
%         plot(t,Gro_Curv(:,j));
%         hold on;
%     end
%     ylim([1E-3 0.5]);
%     subplot(3,7,1);
%     xlabel('Time/min');ylabel('growth crv');
% end

% for i=1:7
%     subplot(3,7,(i+7));
%     for j=(1+(i-1)*12):(i*12)
%         plot(t,RGR(:,j));
%         hold on;
%     end
%     
%     ylim([-2E-5 4E-5]);
%     subplot(3,7,8);
%     xlabel('Time/min');ylabel('growth rt');
% end

%% FFT trans
N=length(RGR(:,1));
Nfft=2^nextpow2(N);
Fs=1/(t(2)-t(1));%sampling rate
f=Fs/2*linspace(0,1,1+Nfft/2);% freq resolution

% test what happens to FFT
for i=1:(7*12)
    FGR(:,i)=fft(RGR(:,i),Nfft)/Nfft;% FFT
    freplot=FGR(:,i);
    FGR_plot(:,i)=2*abs(freplot(1:1+Nfft/2));% extract symmetric part
    FGR_plot1(:,i)=FGR_plot(:,i)/max(FGR_plot(:,i));% normalize the frequency
%     ind(:,i)=find(FGR_plot1(:,i));%record freq point and the amplitude for each 12 replicate
%     figure(2);
%     subplot(3,1,1);
%     stem(abs(FGR(:,i)));
%     subplot(3,1,2);
%     stem(abs(freplot));
%     subplot(3,1,3);
%     stem(FGR_plot1(:,i));
end

%% compare the amplitude for strain B FGR_plot(1:12)
%{
the method to defone colorbar is useful
the log scale
the method to define suptitle
%}

figure;
seqB=FGR_plot(:,1);
n=length(seqB);
cmap=parula(n);
for i= 2:12
    seqBB=FGR_plot(:,i);
    subplot(2,6,i);
    for j=1:n
        loglog(seqB(j),seqBB(j),'o','MarkerFaceColor',cmap(j,:),'MarkerEdgeColor',cmap(j,:));
        hold on;
    end
    ylabel(['rep',num2str(i),'  amplitude'],'Fontsize',14);
    set(gca,'LineWidth',2,'Fontsize',14);
    subplot(2,6,7);
    xlabel('rep1 amplitude','Fontsize',14); 
end
suptitle('Comparison of amplitude');
h=suptitle('Comparison of amplitude in logscale');
set(h,'Fontsize',20);
colorbar


%% plot amplitude comparison separately
seqB= FGR_plot(:,1);
n=length(seqB);
cmap=parula(n);
for i =2:12 
    seqBB=FGR_plot(:,i);
    figure(i);
    for j=1:n
        loglog(seqB(j),seqBB(j),'o','MarkerFaceColor',cmap(j,:),'MarkerEdgeColor',cmap(j,:));
        hold on;
    end
    ylabel(['rep',num2str(i),'  amplitude'],'Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',20);
    xlabel('rep1 amplitude','Fontsize',20); 
    title('Comparison of amplitude','Fontsize',28);
    colorbar
    h=colorbar;
    title(h, 'High Frequency');
    ylabel(h,'Frequency Decreases');
    txt1='Low Frequency';txt2='High Frequency';
    x1=1E-6;y1=1E-6;
    x2=1E-7;y2=1E-7;
    text(x1,y1,txt1,'Fontsize',15);hold on;
    text(x2,y2,txt2,'Fontsize',15);    
end


%% calculate relative standard deviation(RSD) of strain B to H

%{
The way to define RSD, evaluate the statistics
%}

Nam=['B' 'C' 'D' 'E' 'F' 'G' 'H'];
n=length(FGR_plot(:,1));
for j=1:7
    A= FGR_plot(:,(j-1)*12+1:(j-1)*12+12);
    Mean=[];
    Var=[];
    Norm=[];
    for i= 1:n
        Mean(i)=mean(A(i,:));
        Var (i)=var(A(i,:));
        Norm(i)=Var(i)^(1/2)/Mean(i);
    end
    figure(12+j);
    stem(f,Norm,'LineWidth',1.5);
    xlabel('freq/Hz','Fontsize',15);ylabel('Coefficient of Variation','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
    title(sprintf('Coefficient of Variation at Different Frequencies for Strain %s',Nam(j)),'Fontsize',20);
    
end
