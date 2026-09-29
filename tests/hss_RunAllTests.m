function hss_RunAllTests
%RUNALLTESTS Run all the unit tests.


hss_TestCreation;
TestTruncate('hss');

hss_TestOperations;

hss_TestVarious;

hss_TestLyapunov;


end
