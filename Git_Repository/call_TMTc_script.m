% Master script to Run TMTc with sample file


% Define the isolation window for MS2 spectra

%---------------------------------------------- Fusion San Jose
% For 4m/z Iso window
mz_value = [-5:0.001:5];
Fit_Iso_Window =@(x)1.101*sin(0.5564*x+1.591) + -2.168*sin(3.008*x-1.733) + 2.242*sin(2.957*x+-1.729);
Rel_Abundance = Fit_Iso_Window(mz_value);
Rel_Abundance(mz_value < -2.5 | mz_value > 2.5)= 0;
Rel_Abundance(Rel_Abundance<0)=0;
Rel_Abundance = Rel_Abundance./max(Rel_Abundance);
%plot(mz_value,Rel_Abundance,'x')
Array_Iso_Window_4p0mz = [-1*mz_value;Rel_Abundance]';




%For 0.75 m/z isolation window
mz_value = [-5:0.001:5];
Fit_Iso_Window = @(x)  1.803*sin(0.2177*x+ 0.1615) + 0.6801*sin(4.116*x+1.659) +  0.04869*sin(22.93*x+0.8223);
Rel_Abundance = Fit_Iso_Window(mz_value);
Rel_Abundance(mz_value < -.45 | mz_value > 0.6)= 0;
Rel_Abundance(Rel_Abundance<0)=0;
Rel_Abundance = Rel_Abundance./max(Rel_Abundance);
%plot(mz_value,Rel_Abundance,'x')
Array_Iso_Window_0p75mz = [-1*mz_value;Rel_Abundance]';

% For 0.5 m/z isolation window
mz_value = [-5:0.001:5];
Fit_Iso_Window = @(x)  1.081*sin(4.39*x+1.489) + 1.274*sin(18.31*x+1.537) + 1.294*sin(18.83*x-1.631);
Rel_Abundance = Fit_Iso_Window(mz_value);
Rel_Abundance(mz_value < -0.3 | mz_value > 0.45)= 0;
Rel_Abundance(Rel_Abundance<0)=0;
Rel_Abundance = Rel_Abundance./max(Rel_Abundance);
%plot(mz_value,Rel_Abundance,'x')
Array_Iso_Window_0p5mz = [-1*mz_value;Rel_Abundance]';  


%---------------------------------------------- Lumos, Robert 09-13-2015

% For 0.5 m/z isolation window on Lumos
mz_value = [-5:0.001:5];
Fit_Iso_Window = @(x)  1.125*sin(3.579*x+1.432) + 0.4842*sin(16.33*x+1.285) + 0.1908*sin(20.11*x+0.8704) + ... 
                    0.6793*sin(17.53*x+4.291);
Rel_Abundance = Fit_Iso_Window(mz_value);
Rel_Abundance(mz_value < -0.4 | mz_value > 0.55)= 0;
Rel_Abundance(Rel_Abundance<0)=0;
Rel_Abundance = Rel_Abundance./max(Rel_Abundance);
%plot(mz_value,Rel_Abundance,'x')
Array_Iso_Window_0p5mz_Lumos = [mz_value;Rel_Abundance]';  

% For 2 m/z isolation window on Lumos
mz_value = [-5:0.001:5];
       a1 =      0.5934;  
       b1 =      0.4233;  
       c1 =       1.579;  
       a2 =      0.5055;  
       b2 =       1.691; 
       c2 =       1.457;  
       a3 =      0.1363;% (-1.81e+06, 1.81e+06)
       b3 =       5.075;%  (-2.922e+05, 2.922e+05)
       c3 =      -4.337;%  (-3.012e+06, 3.012e+06)
       a4 =     0.09227;%  (-3040, 3040)
       b4 =       6.832;%  (-1.377e+04, 1.378e+04)
       c4 =       1.583;%  (-280.7, 283.9)
       a5 =     -0.2557;%  (-1.828e+06, 1.828e+06)
       b5 =       5.199;%  (-6.818e+05, 6.818e+05)
       c5 =       1.681;%  (-3.102e+05, 3.102e+05)
       a6 =     0.06395;%  (-8.744, 8.872)
       b6 =       10.45;%  (-48.01, 68.9)
       c6 =      -1.707;%  (-8.014, 4.6)
       a7 =      0.1009;%  (-2.395e+04, 2.395e+04)
       b7 =       3.734;%  (-1.549e+05, 1.549e+05)
       c7 =      -1.413;%  (-3.007e+04, 3.007e+04)
       a8 =      0.0105;%  (-18.59, 18.61)
       b8 =       9.415;%  (-240.6, 259.5)
       c8 =      -4.086;%  (-1431, 1423)
Fit_Iso_Window = @(x)   a1*sin(b1*x+c1) + a2*sin(b2*x+c2) + a3*sin(b3*x+c3) + ...
                    a4*sin(b4*x+c4) + a5*sin(b5*x+c5) + a6*sin(b6*x+c6) + ...
                    a7*sin(b7*x+c7) + a8*sin(b8*x+c8);
                
% For 1p4 isolation window on Lumos 

mz_value = [-5:0.001:5];

       a1 =      0.7095 ;
       b1 =      0.8327 ;
       c1 =       1.591 ;
       a2 =       8.143 ;
       b2 =        4.23 ;
       c2 =       1.566 ;
       a3 =      0.5766 ;
       b3 =       10.11 ;
       c3 =      -1.547 ;
       a4 =       7.791 ;
       b4 =       4.295 ;
       c4 =       4.706 ;
       a5 =      0.5433 ;
       b5 =       10.41 ;
       c5 =       1.599 ;
Fit_Iso_Window = @(x)   a1*sin(b1*x+c1) + a2*sin(b2*x+c2) + a3*sin(b3*x+c3) + ...
                    a4*sin(b4*x+c4) + a5*sin(b5*x+c5);
                
Rel_Abundance = Fit_Iso_Window(mz_value);
Rel_Abundance(mz_value < -0.9 | mz_value > 0.9)= 0;
Rel_Abundance(Rel_Abundance<0)=0;
Rel_Abundance = Rel_Abundance./max(Rel_Abundance);
%plot(mz_value,Rel_Abundance,'x')
Array_Iso_Window_1p4mz_Lumos = [mz_value;Rel_Abundance]';  

%Isolation window 0p6 Lumos 

       a1 =      0.4643  ;
       b1 =       1.104  ;
       c1 =       1.791  ;
       a2 =       0.726  ;
       b2 =       3.908 ;
       c2 =       1.493  ;
       a3 =        17.3  ;
       b3 =          15  ;
       c3 =      -1.612  ;
       a4 =      0.7116  ;
       b4 =       13.61  ;
       c4 =        1.47  ;
       a5 =       16.62  ;
       b5 =       15.06  ;
       c5 =       1.531 ;
Fit_Iso_Window = @(x)   a1*sin(b1*x+c1) + a2*sin(b2*x+c2) + a3*sin(b3*x+c3) + ...
                    a4*sin(b4*x+c4) + a5*sin(b5*x+c5);
                
Rel_Abundance = Fit_Iso_Window(mz_value);
Rel_Abundance(mz_value < -0.55 | mz_value > 0.55)= 0;
Rel_Abundance(Rel_Abundance<0)=0;
Rel_Abundance = Rel_Abundance./max(Rel_Abundance);
plot(mz_value,Rel_Abundance,'x')
Array_Iso_Window_0p6mz_Lumos = [mz_value;Rel_Abundance]';        
       

%Call the master program
for filenames = {'v11104_Ynmin1_yeast_only.csv'}
results = demo_TMTc_iso_Window(filenames{1},6,[1,1,1,1,1],1,2000,Array_Iso_Window_4p0mz,0.01);
saving_results_name = strsplit(filenames{1},'.')
save(saving_results_name{1},'results')
clear results
end




