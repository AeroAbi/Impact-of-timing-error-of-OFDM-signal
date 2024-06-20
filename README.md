# Impact-of-timing-error-of-OFDM-signal
For a channel with multipath delay spread, the received signal is a summation of the
transmitted signal with different (complex) gains and delays and there are numerous factors
leading to timing errors in OFDM such as multipath effects, Cyclic prefix imperfections, clock
synchronization errors etc;
Considering 3 cases of FFT misalignment with Cyclic Prefix, the following observations are
recorded
Case 1: The FFT window starts in the range [1, L-1]. This case represents a situation where the
FFT window starts before the cyclic path, resulting in interference from the previous OFDM
symbol. In this case, the received signal will be disturbed due to overlap with previous signal.
Case 2: The FFT window starts in the range [L, D+1]. This case is the ideal situation where
the FFT window starts at the very end of the cyclic path, protecting against interference. In this
case, there should be no interference, and the received signal should not be free of distortion
caused by interference between signals.
Case 3: The FFT window starts between [D+2, D+N]. In this case, the FFT window starts
after the cyclic prefix, adding a portion of the cyclic prefix to the FFT window. This ideally results
in interference from the cyclic path, affecting the received signal. However, compared to Case 1,
the interference should be less, as the circular path is designed to reduce the interference
