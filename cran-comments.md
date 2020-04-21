## Test environments
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (rhub)
* Fedora Linux, R-devel, clang, gfortran (rhub)
* Ubuntu Linux 16.04 LTS, R-release, GCC
* x86_64-w64-mingw32 (64-bit), R-devel (win-builder)
* x86_64-w64-mingw32 (64-bit), R 3.6.3 (win-builder)
* windows-latest, r: 'devel'
* macOS-latest, r: 'devel'
* ubuntu-16.04, r: '3.6'
* xenial, r: "devel"

## R CMD check results

0 errors | 0 warnings | 1 note | 0 skipped

NOTE: Found the following (possibly) invalid DOIs:
  DOI: 10.2202/1557-4679.1098
    From: DESCRIPTION
    Status: libcurl error code 35:
    	schannel: next InitializeSecurityContext failed: SEC_E_ILLEGAL_MESSAGE (0x80090326) - This error usually occurs when a fatal SSL/TLS alert is received (e.g. handshake failed).
    Message: Error
  DOI: 10.2307/2533554
    From: DESCRIPTION
    Status: Forbidden
    Message: 403