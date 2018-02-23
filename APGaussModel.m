fs = 10e3;
t_ap = -1e-3:1/fs:1e-3 - (1/fs);
AP_gauss_template = evalgauss(t_ap,0,6e-8);