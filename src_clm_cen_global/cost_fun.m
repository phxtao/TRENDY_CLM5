function cost = cost_fun(obs_tot_soc, soc_stock_mod)
obs_var = 0.60; % as indicated in WoSIS dataset
soc_layer_weight = 1;
cost = 1*nansum(soc_layer_weight.*((((soc_stock_mod).^(1) - (obs_tot_soc).^(1)).^2)...
    ./(2.*((obs_var*(obs_tot_soc).^(1)).^2))));

end






