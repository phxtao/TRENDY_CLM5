function [cpools_transient, soc_stock_transient, litter_stock_transient, cwd_stock_transient, tot_hr] = fun_transient_trendy_simu(kp, year_start, year_end, initial_cpool, frocing_steady_state, frocing_transient)
%-------------------------------------
% constants
%-------------------------------------
year_size = sum(year_end - year_start)+length(year_start);
time_loop = year_end(1) - year_start(1) + 1;
% time_gap = year_start(2) - year_start(1) + 1;
% if the simulation is default CLM5
is_default = 0;

use_beta = 1;
month_num = 12;

global kelvin_to_celsius
kelvin_to_celsius = 273.15;

global use_vertsoilc npool npool_vr n_soil_layer days_per_year secspday
use_vertsoilc = 1 ; % whether or not use vertical maxing part
npool = 7;  % number of pools if no vertical
npool_vr = 140; % number of pools if vertical
n_soil_layer = 20;  % number of soil layers
days_per_year = 365;
secspday = 24*60*60;
% dt = secspday*30;
days_per_month = [31, 30, 31, 28, 31, 30, 31, 31, 30, 31, 30, 31];

global max_altdepth_cryoturbation max_depth_cryoturb
max_altdepth_cryoturbation = 2;
max_depth_cryoturb = 3;

global dz dz_node zisoi zsoi
% width between two interfaces
dz = [2.000000000000000E-002, 4.000000000000000E-002, 6.000000000000000E-002, ...
    8.000000000000000E-002, 0.120000000000000, 0.160000000000000, ...
    0.200000000000000, 0.240000000000000, 0.280000000000000, ...
    0.320000000000000, 0.360000000000000, 0.400000000000000, ...
    0.440000000000000, 0.540000000000000, 0.640000000000000, ...
    0.740000000000000, 0.840000000000000, 0.940000000000000, ...
    1.04000000000000, 1.14000000000000, 2.39000000000000, ...
    4.67553390593274, 7.63519052838329, 11.1400000000000, ...
    15.1154248593737]';

% depth of the interface
zisoi = [2.000000000000000E-002, 6.000000000000000E-002, ...
    0.120000000000000, 0.200000000000000, 0.320000000000000, ...
    0.480000000000000, 0.680000000000000, 0.920000000000000, ...
    1.20000000000000, 1.52000000000000, 1.88000000000000, ...
    2.28000000000000, 2.72000000000000, 3.26000000000000, ...
    3.90000000000000, 4.64000000000000, 5.48000000000000, ...
    6.42000000000000, 7.46000000000000, 8.60000000000000, ...
    10.9900000000000, 15.6655339059327, 23.3007244343160, ...
    34.4407244343160, 49.5561492936897]';
zisoi_0 = 0;
% depth of the node
zsoi = [1.000000000000000E-002, 4.000000000000000E-002, 9.000000000000000E-002, ...
    0.160000000000000, 0.260000000000000, 0.400000000000000, ...
    0.580000000000000, 0.800000000000000, 1.06000000000000, ...
    1.36000000000000, 1.70000000000000, 2.08000000000000, ...
    2.50000000000000, 2.99000000000000, 3.58000000000000, ...
    4.27000000000000, 5.06000000000000, 5.95000000000000, ...
    6.94000000000000, 8.03000000000000, 9.79500000000000, ...
    13.3277669529664, 19.4831291701244, 28.8707244343160, ...
    41.9984368640029]';

% depth between two node
dz_node = zsoi - [0; zsoi(1:end-1)];

%-------------------------------------
% steady state forcing
%-------------------------------------
% 
% nbedrock_steady_state = frocing_steady_state.nbedrock;
% sand_vector_steady_state = nan;
% 
% npp_mean_steady_state = frocing_steady_state.NPP;
% 
% altmax_current_profile_steady_state = frocing_steady_state.ALTMAX;
% altmax_lastyear_profile_steady_state = frocing_steady_state.ALTMAX_LAST_YEAR;
% 
% xio_steady_state = frocing_steady_state.O_SCALAR;
% xin_steady_state = frocing_steady_state.FPI;
% 
% xit_steady_state = frocing_steady_state.T_SCALAR;
% xiw_steady_state = frocing_steady_state.W_SCALAR;
% 
% npp_annual_steady_state = sum(npp_mean_steady_state.*repmat(days_per_month, [1, length(npp_mean_steady_state)/month_num])*secspday)/(length(npp_mean_steady_state)/month_num); % unit: gc/m2/year
% 
%-------------------------------------
% transient simu forcing
%-------------------------------------

nbedrock_transient = frocing_transient.nbedrock;
sand_vector_transient = nan;

npp_mean_transient = frocing_transient.NPP;

altmax_current_profile_transient = frocing_transient.ALTMAX;
altmax_lastyear_profile_transient = frocing_transient.ALTMAX_LAST_YEAR;

xio_transient = frocing_transient.O_SCALAR;
xin_transient = frocing_transient.FPI;

xit_transient = frocing_transient.T_SCALAR;
xiw_transient = frocing_transient.W_SCALAR;

%-------------------------------------
% define parameters to be optimised
%-------------------------------------
% diffusion (bioturbation) 10^(-4) (m2/yr)
bio = kp(1)*(5*10^(-4) - 3*10^(-5)) + 3*10^(-5);
% cryoturbation 5*10^(-4) (m2/yr)
cryo = kp(2)*(16*10^(-4) - 3*10^(-5)) + 3*10^(-5);
% parameters used in vertical discretization of carbon inputs 10 (metre)
efolding = kp(3)*(1 - 0) + 0;
% turnover time of CWD (yr) 3.3333
tau4cwd = kp(4)*(6 - 1) + 1;
% tau for metabolic litter (yr) 0.0541
tau4l1 = kp(5)*(0.11 - 0) + 0;
% tau for cellulose litter (yr) 0.2041
tau4l2 = kp(6)*(0.3 - 0.1) + 0.1;
% tau for lignin litter (yr)
tau4l3 = tau4l2;
% tau for fast SOC (yr) 0.1370
tau4s1 = kp(7)*(0.5 - 0) + 0;
% tau for slow SOC (yr) 5
tau4s2 = kp(8)*(10 - 1) + 1;
% tau for passive SOC (yr) 222.222
tau4s3 = kp(9)*(400 - 20) + 20;

% fraction from l1 to s2, 0.45
fl1s1 = kp(10)*(0.8 - 0.1) + 0.1;
% fraction from l2 to s1, 0.5
fl2s1 = kp(11)*(0.8 - 0.2) + 0.2;
% fraction from l3 to s2, 0.5
fl3s2 = kp(12)*(0.8 - 0.2) + 0.2;
% fraction from s1 to s2, sand dependeted
fs1s2 = kp(13)*(0.4 - 0) + 0;
% fraction from s1 to s3, sand dependeted
fs1s3 = kp(14)*(0.1 - 0) + 0;
% fraction from s2 to s1, 0.42
fs2s1 = kp(15)*(0.74 - 0.1) + 0.1;
% fraction from s2 to s3, 0.03
fs2s3 = kp(16)*(0.1 - 0) + 0;
% fraction from s3 to s1, 0.45
fs3s1 = kp(17)*(0.9 - 0) + 0;
% fraction from cwd to l2, 0.76
fcwdl2 = kp(18)*(1 - 0.5) + 0.5;
% beta to describe the shape of vertical profile
% beta = 0.95;
beta = kp(19)*(0.9999 - 0.85) + 0.85;
% allocation to different litter pools
p4l1 = kp(20)*(1 - 0) + 0;
p4l2 = kp(21)*(1 - 0) + 0;
p4l3 = kp(22)*(1 - 0) + 0;
p4cwd = 1 - p4l1 - p4l2 - p4l3;

adv = 0; % parameter for advection (m/yr)

% %% steady state simulation
% %----------------------------------------
% % Environmental Scalar (Xi)
% %----------------------------------------
% xit = nan(n_soil_layer, month_num);
% xiw = nan(n_soil_layer, month_num);
% xin = nan(n_soil_layer, month_num);
% xio = nan(n_soil_layer, month_num);
% 
% for imonth = 1:month_num
%     xit(:, imonth) = mean(xit_steady_state(:, imonth:month_num:end), 2);
%     xiw(:, imonth) = mean(xiw_steady_state(:, imonth:month_num:end), 2);
%     xin(:, imonth) = mean(xin_steady_state(:, imonth:month_num:end), 2);
%     xio(:, imonth) = mean(xio_steady_state(:, imonth:month_num:end), 2);
% end
% 
% %----------------------------------------
% % Triangle Matrix, A Matrix and K Matrix
% %----------------------------------------
% sand_vector = mean(sand_vector_steady_state, 2, 'omitnan');
% % Allocation matrix
% a_ma = a_matrix(is_default, fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2, sand_vector);
% 
% kk_ma_middle = nan(npool_vr, npool_vr, month_num);
% tri_ma_middle = nan(npool_vr, npool_vr, month_num);
% for imonth = 1:month_num
%     % decomposition matrix
%     monthly_xit = xit(:, imonth);
%     monthly_xiw = xiw(:, imonth);
%     monthly_xio = xio(:, imonth);
%     monthly_xin = xin(:, imonth);
%     
% 	monthly_xiw(monthly_xiw < 0.05) = 0.05;
% 
% 	kk_ma_middle(:, :, imonth) = kk_matrix(monthly_xit, monthly_xiw, monthly_xio, monthly_xin, efolding, tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3);
%     % tri matrix
%     monthly_nbedrock = nbedrock_steady_state;
%     monthly_altmax_current_profile = mean(altmax_current_profile_steady_state(imonth:month_num:end));
%     monthly_altmax_lastyear_profile = mean(altmax_lastyear_profile_steady_state(imonth:month_num:end));
%     tri_ma_middle(:, :, imonth) = tri_matrix(monthly_nbedrock, monthly_altmax_current_profile, monthly_altmax_lastyear_profile, bio, adv, cryo);
% end
% tri_ma = mean(tri_ma_middle, 3, 'omitnan');
% kk_ma = mean(kk_ma_middle, 3, 'omitnan');
% 
% %----------------------------------------
% % Vertical Profile
% %----------------------------------------
% % in the original beta model in Jackson et al 1996, the unit for the depth
% % of the soil is cm (dmax*100)
% 
% m_to_cm = 100;
% vertical_prof = nan(n_soil_layer, 1);
% 
% if mean(altmax_lastyear_profile_steady_state) > 0
%     for j = 1:n_soil_layer
%         if j == 1
%             vertical_prof(j) = (beta^((zisoi_0)*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
%         else
%             vertical_prof(j) = (beta^((zisoi(j - 1))*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
%         end
%     end
% else
%     vertical_prof(1) = 1/dz(1);
%     vertical_prof(2:end) = 0;
% end
% 
% vertical_input = dz(1:n_soil_layer).*vertical_prof/sum(vertical_prof.*dz(1:n_soil_layer));
% %----------------------------------------
% % Analytical Solution of SOC
% %----------------------------------------
% matrix_in=nan(140,1);
% % total input amount
% input_tot_cwd = p4cwd*npp_annual_steady_state/days_per_year; % (gc/m2/day)
% input_tot_litter1 = p4l1*npp_annual_steady_state/days_per_year; % (gc/m2/day)
% input_tot_litter2 = p4l2*npp_annual_steady_state/days_per_year; % (gc/m2/day)
% input_tot_litter3 = p4l3*npp_annual_steady_state/days_per_year; % (gc/m2/day)
% 
% % redistribution by beta
% matrix_in(1:20,1) = input_tot_cwd*vertical_input./dz(1:n_soil_layer); % litter input gc/m3/day
% matrix_in(21:40,1) = input_tot_litter1*vertical_input./dz(1:n_soil_layer);
% matrix_in(41:60,1) = input_tot_litter2*vertical_input./dz(1:n_soil_layer);
% matrix_in(61:80,1) = input_tot_litter3*vertical_input./dz(1:n_soil_layer);
% matrix_in(81:140,1) = 0;
% 
% input_matrix(1:20,1) = p4cwd.*vertical_input./dz(1:n_soil_layer);
% input_matrix(21:40,1) = p4l1.*vertical_input./dz(1:n_soil_layer);
% input_matrix(41:60,1) = p4l2.*vertical_input./dz(1:n_soil_layer);
% input_matrix(61:80,1) = p4l3.*vertical_input./dz(1:n_soil_layer);
% input_matrix(81:140,1) = 0;
% 
% lastwarn('');
% 
% cpool_steady_state = (a_ma*kk_ma-tri_ma)\(-matrix_in);
% cpool_steady_state(cpool_steady_state < 0) = 0;
% [~, msgid] = lastwarn;
% if strcmp(msgid, 'MATLAB:singularMatrix') == 1 || strcmp(msgid, 'MATLAB:nearlySingularMatrix') == 1
%     cpools_transient = nan;
%     soc_stock_transient = nan;
%     litter_stock_transient = nan;
%     cwd_stock_transient = nan;
%     tot_hr = nan;
%     
%     return
% end
% 
% initial_cpool(:, 1) = cpool_steady_state;
%% transient simulation
% soc related variables
cpools_transient = nan(year_size*month_num, npool*n_soil_layer);
soc_stock_transient = nan(year_size*month_num, 1);
litter_stock_transient = nan(year_size*month_num, 1);
cwd_stock_transient = nan(year_size*month_num, 1);
tot_hr = nan(year_size*month_num, 1);
% env. scalars
xit = nan(n_soil_layer, 1);
xiw = nan(n_soil_layer, 1);
xio = nan(n_soil_layer, 1);
xin = nan(n_soil_layer, 1);

for iyear = 1 : year_size
    for imonth = 1 : month_num
        if imonth == 1 && mod(iyear-1, time_loop) == 0
            current_cpool = initial_cpool(:, round((iyear-1)/time_loop)+1);
        end
        
        month_loc = (iyear-1)*month_num+imonth;
        
        xit(:, 1) = xit_transient(:, month_loc);
        xiw(:, 1) = xiw_transient(:, month_loc);
        xin(:, 1) = xin_transient(:, month_loc);
        xio(:, 1) = xio_transient(:, month_loc);
        
        %----------------------------------------
        % Triangle Matrix, A Matrix and K Matrix
        %----------------------------------------
        % Allocation matrix
        a_ma = a_matrix(is_default, fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2, sand_vector_transient);
        
        kk_ma = kk_matrix(xit, xiw, xio, xin, efolding, tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3);
        % tri matrix
        tri_ma = tri_matrix(nbedrock_transient, altmax_current_profile_transient(month_loc), altmax_lastyear_profile_transient(month_loc), bio, adv, cryo);
        
        %----------------------------------------
        % Vertical Profile
        %----------------------------------------
        % in the original beta model in Jackson et al 1996, the unit for the depth
        % of the soil is cm (dmax*100)
        if is_default == 0
            m_to_cm = 100;
            vertical_prof = nan(n_soil_layer, 1);
            
            altmax_lastyear_profile = altmax_lastyear_profile_transient(month_loc);
            if altmax_lastyear_profile > 0
                for j = 1:n_soil_layer
                    if j == 1
                        vertical_prof(j) = (beta^((zisoi_0)*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
                    else
                        vertical_prof(j) = (beta^((zisoi(j - 1))*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
                    end
                end
            else
                vertical_prof(1) = 1/dz(1);
                vertical_prof(2:end) = 0;
            end
            
            vertical_input = dz(1:n_soil_layer).*vertical_prof/sum(vertical_prof.*dz(1:n_soil_layer));
        end
        %----------------------------------------
        % input vector
        %----------------------------------------
        
        matrix_in=nan(140,1);
        % total input amount
        input_tot_cwd = p4cwd*npp_mean_transient(month_loc)*secspday; % (gc/m2/day)
        input_tot_litter1 = p4l1*npp_mean_transient(month_loc)*secspday; % (gc/m2/day)
        input_tot_litter2 = p4l2*npp_mean_transient(month_loc)*secspday; % (gc/m2/day)
        input_tot_litter3 = p4l3*npp_mean_transient(month_loc)*secspday; % (gc/m2/day)
        carbon_input = npp_mean_transient(month_loc)*secspday; % (gc/m2/day)
        
        % redistribution by beta
        matrix_in(1:20,1) = input_tot_cwd*vertical_input./dz(1:n_soil_layer); % litter input gc/m3/day
        matrix_in(21:40,1) = input_tot_litter1*vertical_input./dz(1:n_soil_layer);
        matrix_in(41:60,1) = input_tot_litter2*vertical_input./dz(1:n_soil_layer);
        matrix_in(61:80,1) = input_tot_litter3*vertical_input./dz(1:n_soil_layer);
        matrix_in(81:140,1) = 0;
        
        input_matrix(1:20,1) = p4cwd.*vertical_input./dz(1:n_soil_layer);
        input_matrix(21:40,1) = p4l1.*vertical_input./dz(1:n_soil_layer);
        input_matrix(41:60,1) = p4l2.*vertical_input./dz(1:n_soil_layer);
        input_matrix(61:80,1) = p4l3.*vertical_input./dz(1:n_soil_layer);
        input_matrix(81:140,1) = 0;
        
        %----------------------------------------
        % soc stepwise iteration
        %----------------------------------------
        cpools_next = (matrix_in + (a_ma*kk_ma-tri_ma)*current_cpool)*days_per_month(imonth) + current_cpool;
        
        %         residence_time = ((a_ma*kk_ma - tri_ma)\(-input_matrix))/days_per_year;
        %         residence_time_layer =  sum([residence_time(81:100), residence_time(101:120), residence_time(121:140)], 2); % unit yr
        %         residence_time_30cm = sum(residence_time_layer(1:4).*dz(1:4)) + residence_time_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        %         residence_time_100cm = sum(residence_time_layer(1:8).*dz(1:8)) + residence_time_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        %         residence_time_200cm = sum(residence_time_layer(1:11).*dz(1:11)) + residence_time_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        %         residence_time_total = sum(residence_time_layer.*dz(1:n_soil_layer));
        
        soc_layer = sum([cpools_next(81:100), cpools_next(101:120), cpools_next(121:140)], 2); % unit: g c/m3
        %         soc_stock_30cm = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4)); % unit gc/m2
        %         soc_stock_100cm = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8)); % unit gc/m2
        %         soc_stock_200cm = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11)); % unit gc/m2
        soc_stock_total = sum(soc_layer.*dz(1:n_soil_layer)); % unit gc/m2
        
        litter_layer = sum([cpools_next(21:40), cpools_next(41:60), cpools_next(61:80)], 2); % unit: g c/m3
        cwd_layer = cpools_next(1:20); % unit: g c/m3
        
        litter_stock_total = sum(litter_layer.*dz(1:n_soil_layer)); % unit gc/m2
        cwd_stock_total = sum(cwd_layer.*dz(1:n_soil_layer)); % unit gc/m2
        %----------------------------------------
        % respired soc
        %----------------------------------------
        r_vector = -sum(a_ma*kk_ma, 1); % unit: day-1
        
        % respired soc
        respired_cpools = r_vector'.*(current_cpool)*days_per_month(imonth); % unit: gc m-3 month-1
        respired_layers = [respired_cpools(1:20), respired_cpools(21:40), respired_cpools(41:60), respired_cpools(61:80), ...
            respired_cpools(81:100), respired_cpools(101:120), respired_cpools(121:140)];
        respired_layers = sum(respired_layers, 2); % unit: gc m-3 month-1
        respired_total = sum(respired_layers.*dz(1:n_soil_layer))/days_per_month(imonth)/secspday; % unit: gc m-2 sec-1
        
        %----------------------------------------
        % carbon pool update
        %----------------------------------------
        cpools_transient((iyear-1)*12+imonth, :) = cpools_next;
        soc_stock_transient((iyear-1)*12+imonth, :) = soc_stock_total;
        
        litter_stock_transient((iyear-1)*12+imonth, :) = litter_stock_total;
        cwd_stock_transient((iyear-1)*12+imonth, :) = cwd_stock_total;
        tot_hr((iyear-1)*12+imonth, :) = respired_total;
        
        current_cpool = cpools_next;
        
    end % imonth = ...
end % iyear = ...

end  % function

% ----- CENTURY T response function
function catanf_results = catanf(t1)
catanf_results = 11.75 +(29.7 / pi) * atan( pi * 0.031  * ( t1 - 15.4 ));
end


