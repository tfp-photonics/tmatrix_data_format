function [epsilon] = epsilon_Si(wavelength)
%% epsilon_Si
% Dielectric function of silicon in the UV-vis-nearIR region from tabulated
% values (D. E. Aspnes and A. A. Studna. Phys. Rev. B 27, 985-1009 (1983))
% 10.1103/PhysRevB.27.985


d = [206.6000    1.0100    2.9090;
    210.1000    1.0830    2.9820;
    213.8000    1.1330    3.0450;
    217.5000    1.1860    3.1200;
    221.4000    1.2470    3.2060;
    225.4000    1.3400    3.3020;
    229.6000    1.4710    3.3660;
    233.9000    1.5790    3.3530;
    238.4000    1.5890    3.3540;
    243.1000    1.5710    3.4290;
    248.0000    1.5700    3.5650;
    253.0000    1.5970    3.7490;
    258.3000    1.6580    3.9790;
    263.8000    1.7640    4.2780;
    269.5000    1.9880    4.6780;
    275.5000    2.4520    5.0820;
    281.8000    3.1200    5.3440;
    288.3000    4.0870    5.3950;
    295.2000    4.8880    4.6390;
    302.4000    5.0200    3.9790;
    310.0000    5.0100    3.5860;
    317.9000    5.0160    3.3460;
    326.3000    5.0650    3.1820;
    335.1000    5.1560    3.0580;
    344.4000    5.2960    2.9870;
    354.2000    5.6100    3.0140;
    364.7000    6.5220    2.7050;
    375.7000    6.7090    1.3200;
    387.5000    6.0620    0.6300;
    399.9000    5.5700    0.3870;
    413.3000    5.2220    0.2690;
    427.5000    4.9610    0.2030;
    442.8000    4.7530    0.1630;
    459.2000    4.5830    0.1300;
    476.9000    4.4420    0.0900;
    495.9000    4.3200    0.0730;
    516.6000    4.2150    0.0600;
    539.1000    4.1230    0.0480;
    563.6000    4.0420    0.0320;
    590.4000    3.9690    0.0300;
    619.9000    3.9060    0.0220;
    652.5000    3.8470    0.0160;
    688.8000    3.7960    0.0130;
    729.3000    3.7520    0.0100;
    774.9000    3.7140    0.0080;
    826.6000    3.6730    0.0050;
    1200        3.6730    0.0050]; % added for extended interpolation

epsr  = interp1(d(:,1), d(:,2), wavelength, "linear");
epsi  = interp1(d(:,1), d(:,3), wavelength, "linear");

epsilon = epsr + 1i*epsi;

end
