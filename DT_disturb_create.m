function Fdisturbance_mm = DT_disturb_create(p, nsample)
% load output additive disturbances
% sampling frequency: 10 kHz
% duration: 5 s -- 50K samples
% disturb magnitude ; micron
%
% disturb: spread spectrum around 20--30Hz, 
% then a few harmonics of the 50 Hz electric power feed-in
% 
    load Y_Add_Disturb.mat Fdisturbance 
    
    Fdisturbance = Fdisturbance(1:p,1:nsample); 
    
    % Fdisturbance in micron
    Fdisturbance_mm = Fdisturbance * 1e-3; 
    % conversion from micron [mu-m] to millimeters [mm]

end