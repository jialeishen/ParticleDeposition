function [vd] = vd_empirical(dp, theta, u0)
% Empirical equation for particle dry deposition on surface
% Source: You, R., Zhao, B., Chen, C., 2012. Developing an empirical equation for modeling particle deposition velocity onto inclined surfaces in indoor environments. Aerosol Sci. Technol. 46, 1090–1099. 
% https://doi.org/10.1080/02786826.2012.695096
% vd: deposition velocity [m/s]
% dp: particle diameter [um]
% theta: inclination angle [deg]
% u0: friction velocity [cm/s]

for i = 1:length(dp)
    if dp(i) < (0.0512.*u0.^0.4227)
        vd(i) = (5.15e-8.*u0 - 5.63e-11).*dp(i).^(-1.263);
    elseif (dp(i) > (0.3577.*cosd(theta).^(-0.41))) && (cosd(theta) >= 0)
        vd(i) = 3.7e-5.*dp(i).^1.9143.*cosd(theta);
    elseif (log10(dp(i)) > -0.941 + 0.796.*log10(u0) + 0.333.*cosd(theta) + 0.184.*log10(u0).*cosd(theta) - 0.011.*(log10(u0)).^2 + 0.15.*(cosd(theta)).^2) && (cosd(theta) < 0)
        vd(i) = 0;
    else
        f = -6.026 + 0.116.*log10(u0) + 1.837.*cosd(theta) + 1.079.*log10(dp(i)) - 0.653.*log10(u0).*cosd(theta) - 0.89.*log10(u0).*log10(dp(i)) + 1.484.*log10(dp(i)).*cosd(theta) + 0.17.*(log10(u0)).^2 - 0.074.*(cosd(theta)).^2 + 1.076.*(log10(dp(i))).^2;
        vd(i) = 10.^f;
    end
end

end

