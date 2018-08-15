# ElaCov
Matlab version of code is recommended.

For R, the myMatrix2tmpCppCoreForR.dll is genenrated on a windows-platform PC. The user
needs to regenerate the .dll file if the code doesn't work.

Use the following command in terminal to generate DLL for windows or shared library for Mac and Linux

R CMD SHLIB myMatrix2tmpCppCoreForR.c

Then for Mac or Linux user, you also need to change 'myMatrix2tmpCppCoreForR.dll' in line 113 and 124 in ElaCov.R to 'myMatrix2tmpCppCoreForR.so', after using 'R CMD SHLIB'.
