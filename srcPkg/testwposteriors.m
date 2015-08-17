randn('state', 12345);

rate = 2048;

channelNames = { 'H1:foo', 'L1:foo', 'V1:foo' };

x = {};
xSw = {};

t = ((1:rate)-.5)/rate;
envelope = 6*exp(-((t - 0.5).^2)/(2 * 0.001.^2));
phase = 2 * pi * (t - 0.5) * 256;
wc = envelope .* cos(phase);
ws = envelope .* sin(phase);
wSw = repmat(wc * wc', 1, length(channelNames));

figure(1);
hold off;

marker = { 'r-', 'g-', 'm-' };

taus = [];

for k = 1:length(channelNames)
    x{k} = randn(1, rate);
    
    if (k)
        x{k} = x{k} + [ wc((rate/2+1):rate) , wc(1:(rate/2)) ];
    end
        
    xSw{k} = real(ifft(fft(x{k}) .* conj(fft(wc)))) ...
       + i * real(ifft(fft(x{k}) .* conj(fft(ws))));
   
   plot(t, abs(xSw{k}), marker{k});
   hold on;

   taus = [taus, 0.5 - 0.2.^k, 0.5+0.2.^k];

end

taus = [taus, 0, 1, 0.0001];


figure(2);

% mean(real(xSw) .* real(xSw)) / wSw

[theta, phi] = wsinusoidalprojection(180);
directions = [ theta(:)' ; phi(:)' ];

%taus = rate/2 + ((-127):128);

[logSkymap] = wposteriors(channelNames, rate, wSw, xSw, directions, taus);

whybridskyplot([theta(:) , phi(:) , logSkymap(:)]);

logSignal = log(mean(exp(logSkymap - max(logSkymap)))) + max(logSkymap)

