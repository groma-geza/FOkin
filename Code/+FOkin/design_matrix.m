function A = design_matrix(time, t0, tau, fwhm)
% Design matrix for exponential decay.

m = length(time);
n = length(tau);
if fwhm == 0
    [taumatr,time] = meshgrid(tau,time);
    A=exp(-(time-t0)./taumatr);
    A(time < t0)= 0;
else
    A = zeros(m, n);
    sigma = fwhm/2/sqrt(2*log(2));
    for ii=1:n
        col = conv_gauss_exp(time, t0, sigma, tau(ii));
        A(:,ii) = col';
    end
end
end

%%
function y = conv_gauss_exp(t, t0, sigma, tau)

% Calculates the convolution of a Gaussian and an exponential
% Basic equations in LaTeX format
%
% \[\begin{array}{*{20}{l}}
% {\begin{array}{*{20}{c}}
% {gauss(\sigma ,t) = \frac{1}{{\sqrt {2\pi } \sigma }}{e^{ - \frac{{{t^2}}}{{2{\sigma ^2}}}}},\quad \quad \sigma  > 0}\\
% {}\\
% {fwh{m_{gauss(\sigma ,t)}} = \sigma 2\sqrt {2\log (2)} }\\
% {}\\
% {gauss(\sigma ,t - {t_0}) \otimes 1(t){e^{ - \frac{t}{\tau }}} = \frac{1}{2}{e^{\frac{1}{{2\tau }}\left[ { - 2\left( {t - {t_0}} \right) + \frac{{{\sigma ^2}}}{\tau }} \right]}}\left[ {1 + erf\left( {\frac{{t - {t_0} - \frac{{{\sigma ^2}}}{\tau }}}{{\sqrt 2 \sigma }}} \right)} \right]}
% \end{array}}\\
% {1 + erf(x) \equiv erfc( - x)}
% \end{array}\]


% Using erfc(-x) instead of 1+erf(x) provides more exact results.

y = 0.5*exp(sigma^2/2/tau^2)*exp(-(t-t0)/tau).*erfc(-(t-t0-sigma^2/tau)/sigma/sqrt(2));
end