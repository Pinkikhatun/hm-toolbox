function hodlr_RunAllTests
%hodlr_RUNALLTESTS Run all the unit tests for @hodlr

hodlr_TestCreation;
TestTruncate('hodlr');
hodlr_TestOperations;
hodlr_TestLyapunov;

end

