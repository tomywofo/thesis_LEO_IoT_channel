function time = find_20(theta_max, height, inclination, vis_window)
upper = 0;
lower = -vis_window/2;
test_val = (upper+lower)/2;
estim = elevation_angle(theta_max, height, inclination, test_val);
dif = estim-20*pi/180;
count = 0;
acc = 1e-5;
while abs(dif) > acc
    if dif > 0
        upper = (upper+lower)/2;
    else
        lower = (upper+lower)/2;
    end
    test_val = (upper+lower)/2;
    estim = elevation_angle(theta_max, height, inclination, test_val);
    dif = estim-20*pi/180;
    count = count + 1;
end
time = test_val;