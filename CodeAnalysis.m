subplot(2,1,1)
hold off
Cavity=2;
range=520:620;
Concentration=[369 865 1440 958 225]+0;
zeroindex=30;
alpha1=(CavityCounts{Cavity}(zeroindex,range)-CavityCounts{Cavity}(10,range))./CavityCounts{Cavity}(10,range);
alpha2=(CavityCounts{Cavity}(zeroindex,range)-CavityCounts{Cavity}(13,range))./CavityCounts{Cavity}(13,range);
alpha3=(CavityCounts{Cavity}(zeroindex,range)-CavityCounts{Cavity}(16,range))./CavityCounts{Cavity}(16,range);
alpha4=(CavityCounts{Cavity}(zeroindex,range)-CavityCounts{Cavity}(31,range))./CavityCounts{Cavity}(31,range);
alpha5=(CavityCounts{Cavity}(zeroindex,range)-CavityCounts{Cavity}(34,range))./CavityCounts{Cavity}(34,range);

plot(WaveLength(range), alpha1/Concentration(1),'o')
hold all
plot(WaveLength(range), alpha2/Concentration(2),'o')
%hold all
plot(WaveLength(range), alpha3/Concentration(3),'o')
plot(WaveLength(range), alpha4/Concentration(4),'o')
plot(WaveLength(range), alpha5/Concentration(5),'o')
xlabel('Wavelength (nm)')

subplot(2,1,2)
hold off
plot(WaveLength(range), CavityCounts{Cavity}(9,range))
hold all
plot(WaveLength(range), CavityCounts{Cavity}(10,range))
plot(WaveLength(range), CavityCounts{Cavity}(13,range))
plot(WaveLength(range), CavityCounts{Cavity}(16,range))
plot(WaveLength(range), CavityCounts{Cavity}(31,range))
plot(WaveLength(range), CavityCounts{Cavity}(34,range))
