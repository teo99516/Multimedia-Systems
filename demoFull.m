cd level1
SNR1 = demoAAC1('../data/LicorDeCalandraca.wav', '../data/decodedLicorDeCalandraca2.wav');
fprintf('SNR1 = %.4f dB\n',SNR1);
cd ..
cd level2
SNR2 = demoAAC2('../data/LicorDeCalandraca.wav', '../data/decodedLicorDeCalandraca2.wav');
fprintf('SNR2 = %.4f dB\n',SNR2);
cd ..
cd level3
SNR3 = demoAAC3('../data/LicorDeCalandraca.wav', '../data/decodedLicorDeCalandraca3.wav','../data/AACseq3.mat');
fprintf('SNR3 = %.4f dB\n',SNR3);
cd ..