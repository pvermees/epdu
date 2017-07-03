clear all;
%% load the data
% load white_heliumages;
% data = white_heliumages(:,1)';
% errors = white_heliumages(:,2)';
load bimodal;
% bimodaldata = [10+rand(1,10)*10 70+rand(1,10)*5];
% bimodalerrors = 1+rand(1,20)*5;
data = bimodaldata;
errors = bimodalerrors;

% epdu(data,errors,2,1,1,1,0); % 1 iteration
% exportfig(gcf,'mb=2_ns=1_tb=1_ni=1','format','jpeg','height',2,'Resolution',300);
% close all;
% epdu(data,errors,1,1,1,100,0); % 100 iterations
% exportfig(gcf,'mb=1_ns=1_tb=1_ni=100','format','jpeg','height',2,'Resolution',300);
% close all;
epdu(data,errors,1,100,10,100,0); % 10000 iterations
% exportfig(gcf,'mb=1_ns=1_tb=1_ni=10000','format','jpeg','height',2,'Resolution',300);
% close all;
% epdu(data,errors,0,100,1,1,0); % 100 iterations
% exportfig(gcf,'mb=0_ns=100_tb=1_ni=1','format','jpeg','height',2,'Resolution',300);
% close all;
% epdu(data,errors,0,10000,1,1,0); % 10000 iterations
% exportfig(gcf,'mb=0_ns=10000_tb=1_ni=1','format','jpeg','height',2,'Resolution',300);
% close all;
% epdu(data,errors,0,10000,10,1,0); % 10000 iterations
% exportfig(gcf,'mb=0_ns=10000_tb=10_ni=1','format','jpeg','height',2,'Resolution',300);
% close all;
% epdu(data,errors,0,1,10,10000,0); % 10000 iterations
% exportfig(gcf,'mb=0_ns=1_tb=10_ni=10000','format','jpeg','height',2,'Resolution',300);
% close all;
% epdu(data,errors,0,100,1,100,0); % 10000 iterations
% exportfig(gcf,'mb=0_ns=100_tb=1_ni=100','format','jpeg','height',2,'Resolution',300);
% close all;
% epdu(data,errors,0,100,10,100,0); % 10000 iterations
% exportfig(gcf,'mb=0_ns=100_tb=10_ni=100','format','jpeg','height',2,'Resolution',300);
% 
% epdu(data,errors,0,1,4,1,1); % 4 iterations
% exportfig(gcf,'mb=0_ns=1_tb=4_ni=1_coupled','format','jpeg','height',2,'Resolution',300);
% close all;
% epdu(data,errors,0,1,10,1,1); % 10 iterations
% exportfig(gcf,'mb=0_ns=1_tb=10_ni=1_coupled','format','jpeg','height',2,'Resolution',300);
% close all;
% epdu(data,errors,0,1000,10,1,1); % 10000 iterations
% exportfig(gcf,'mb=0_ns=1000_tb=10_ni=1_coupled','format','jpeg','height',2,'Resolution',300);
% close all;
% epdu(data,errors,0,1,10,1000,1); % 10000 iterations
% exportfig(gcf,'mb=0_ns=1_tb=10_ni=1000_coupled','format','jpeg','height',2,'Resolution',300);
% close all;
% epdu(data,errors,0,50,4,50,1); % 10000 iterations
% exportfig(gcf,'mb=0_ns=50_tb=4_ni=50_coupled','format','jpeg','height',2,'Resolution',300);
% close all;