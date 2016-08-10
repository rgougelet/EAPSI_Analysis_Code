// FileResize.c
// Resize an existing file
//
// [Status, Msg] = FileResize(File, Bytes, [Mode])
// INPUT:
//   File:  Name of an existing file as string, with or without absolute or
//          relative path. Unicode and UNC paths are accepted.
//   Bytes: Number of bytes as numerical scalar.
//          Use UINT64, if the file size exceeds 9.01 PetaBytes.  :-)
//   Mode:  String to describe the mode:
//          'set':  Set the file size to the absolute size.
//          'move': Append or crop bytes relative to the end. If [Bytes] is
//                  negative, the file size is reduced.
//          Optional, default: 'set'.
//
// OUTPUT:
//   Status: Scalar DOUBLE. Optional, if the output is not caught by the caller,
//           this function stops with an error on problems.
//            0: Success
//           -1: Path is a folder or file is write protected
//           -2: File not found
//           -3: No space on device
//           -4: File length is negative
//           -5: File is write protected or in use
//           -6: Unknown proble
//   Msg: String, empty on success, some information in case of problems.
//
// NOTES:
// - Appended bytes are filled by zeros.
// - Resizing a file is more efficient than appending zeros manually.
//
// EXAMPLES:
//   File = tempname; fid = fopen(File, 'w'); fclose(fid);
//   FileResize(File,  8);           % 8 zeros
//   FileResize(File,  8, 'move');   % 16 zeros
//   FileResize(File, -6, 'move');   % 10 zeros
//   delete(File);
//
// COMPILE:
// The C-file must be compiled before it can be used.
// Call "mex -setup" on demand.
//   mex -O FileResize.c
// Linux: consider C99 comments:
//   mex -O CFLAGS="\$CFLAGS -std=c99" FileResize.c
// This function cannot be compiled with LCC2.4 of Matlab 6.5, but with LCC of
// Matlab 2009a/32.
// Pre-compiled Mex: http://www.n-simon.de/mex
// Run the unit-test function uTest_FileResize after compiling.
//
// Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
//         Compiler: LCC3.8, BCC5.5, OWC1.8, MSVC2008/2010
// Assumed Compatibility: higher Matlab versions, Linux and MacOS.
// Author: Jan Simon, Heidelberg, (C) 2012 matlab.THISYEAR(a)nMINUSsimon.de

/*
% $JRev: R-m V:012 Sum:6FBY7sH/tLCj Date:01-Jul-2012 13:39:36 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_FileResize $
% $File: Tools\Mex\Source\FileResize.c $
% History:
% 001: 26-Mar-2012 00:38, First version.
*/

#if defined(__WINDOWS__) || defined(_WIN32) || defined(_WIN64)  // Windows: ----
#  include <windows.h>
#  include <wchar.h>

#  define USE_WCHAR 1

// Under Windows the API and the CRT methods have an equivalent speed, so it
// might be a question of taste:
#  define WINDOWS_API_METHOD 1
#  define WINDOWS_CRT_METHOD 0

#  define POSIX_METHOD       0  // DO NOT ENABLE THIS UNDER WINDOWS

// BCC5.5 forgot the attribute INVALID?! A strange bug.
#  ifndef INVALID_FILE_ATTRIBUTES
#    define INVALID_FILE_ATTRIBUTES ((DWORD) -1)
#  endif

#else  // Linux or MacOS: ------------------------------------------------------
#  define USE_WCHAR 0           // UTF-8 instead of 4-byte CHAR under unix/linux

#  define WINDOWS_API_METHOD 0  // DO NOT ENABLE THIS UNDER LINUX/MacOS
#  define WINDOWS_CRT_METHOD 0  // DO NOT ENABLE THIS UNDER LINUX/MacOS
#  define POSIX_METHOD       1

// Convert Matlab's CHARs, which are UINT16, to UTF-8 unicode strings:
// mxArrayToString() is documented, but the docs forgat to explain the UTF-8
// compatibility. mxArrayToString_UTF8() is not documented, but it works:
// Prototype not needed in 6.5, 2009a, 2011b:
// char *mxArrayToString_UTF8(const mxArray *pa);
#  define MX_ARRAY_TO_STRING mxArrayToString_UTF8
// #  define MX_ARRAY_TO_STRING mxArrayToString

#  include "io64.h"
#endif

// For all operating systems: --------------------------------------------------
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#include <stdio.h>
#include <errno.h>
#include "mex.h"

// Character type for this machine:
#if USE_WCHAR
#define XCHAR wchar_t
#else
#define XCHAR char
#endif

// Assume 32 bit addressing for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T           // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// Error messages do not contain the function name in Matlab 6.5! This is not
// necessary in Matlab 7, but it does not bother:
#define ERR_ID   "JSimon:FileResize:"
#define ERR_HEAD "*** FileResize[mex]: "

// Enumerator for modes:
typedef enum {MODE_BEGIN, MODE_END} Mode_t;

// Error messages;
static const char
        *MSG_bad_file  = "Given path is a directory, or file is read-only.",
        *MSG_no_file   = "File not found.",
        *MSG_disk_full = "No space on device.",
        *MSG_neg_size  = "File length must be positive.",
        *MSG_locked    = "File is write protected or in use.",
        *MSG_unknown   = "Unknown problem.";
static enum {
        ID_ok        =  0,
        ID_bad_file  = -1,
        ID_no_file   = -2,
        ID_disk_full = -3,
        ID_neg_size  = -4,
        ID_locked    = -5,
        ID_unknown   = -6};

// Prototypes:
XCHAR *mxCharToXChar(const mxArray *S);
void Core(XCHAR *FileName, int64_T Distance, Mode_t Mode, int *Status,
          const char **ErrMsg);

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  XCHAR   *FileName;
  int64_T Distance;
  Mode_t  Mode = MODE_BEGIN;
  int     Status = ID_ok;
  const char *ErrMsg = "unknown problem";
            
  // Check number and type of arguments: ---------------------------------------
  if (nrhs != 2 && nrhs != 3) {
     mexErrMsgIdAndTxt(ERR_ID   "BadNInput",
                       ERR_HEAD "2 or 3 inputs required.");
  }
  if (nlhs > 2) {
     mexErrMsgIdAndTxt(ERR_ID   "BadNOutput",
                       ERR_HEAD "2 outputs allowed.");
  }
  
  // Type of input arguments:
  if (!mxIsChar(prhs[0]) || mxIsEmpty(prhs[0])) {
     mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput1",
                       ERR_HEAD "File name must be a non-empty string.");
  }
  if (!mxIsNumeric(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1) {
     mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput2",
                       ERR_HEAD "[Bytes] must be a numerical scalar.");
  }
  
  // Parse 3rd input - check first character only:
  if (nrhs == 3) {
     if (!mxIsChar(prhs[2])) {
        mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput3",
                          ERR_HEAD "3rd input [Mode] must be a string.");
     }
     
     if (!mxIsEmpty(prhs[2])) {
        switch (*(mxChar *) mxGetData(prhs[2])) {
           case L's':
           case L'S':  Mode = MODE_BEGIN;  break;
           case L'm':
           case L'M':  Mode = MODE_END;    break;
           default:
              mexErrMsgIdAndTxt(ERR_ID "UnknownStringInput3",
                              ERR_HEAD "Input [Mode] must be 'set' or 'move'.");
        }
     }
  }
  
  // Obtain file name as wchar:
  FileName = mxCharToXChar(prhs[0]);

  // Get signed distance:
  Distance = (int64_T) mxGetScalar(prhs[1]);
          
  // Call the core function: ---------------------------------------------------
  Core(FileName, Distance, Mode, &Status, &ErrMsg);
  
  // Catch problems: -----------------------------------------------------------
  if (nlhs > 0) {
     plhs[0] = mxCreateDoubleScalar((double) Status);
     if (nlhs > 1) {
        plhs[1] = mxCreateString(ErrMsg);
     }
     
  } else if (Status != 0) {  // Error occurred:
     mexErrMsgIdAndTxt(ERR_ID   "Failed",
                       ERR_HEAD "Cannot resize file: %s\n    File: %ws",
                       ErrMsg, FileName);
  }
  
  // Release memory:
  mxFree(FileName);
    
  return;
}

// =============================================================================
XCHAR *mxCharToXChar(const mxArray *S)
{
  // Convert Matlab string to UTF-16 wchar string on Windows and UTF-8 string on
  // Linux and MacOS.
   
  mwSize Len;
  XCHAR  *W;
  
  // Alternative:
  // if (mbstowcs(W, mxGetData(S), Len) != Len + 1) {
  //   mexErrMsgIdAndTxt(ERR_ID   "BadName",
  //                     ERR_HEAD "Cannot convert name to WCHAR.");
  // }

#if USE_WCHAR  // 2 byte wchar_t for Windows:
  // Reserver wchar string:
  Len = mxGetNumberOfElements(S);
  W   = (XCHAR *) mxMalloc((Len + 1) * sizeof(XCHAR));
  if (W == NULL) {
     mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
                       ERR_HEAD "Cannot get memory for unicode string.");
  }
  
  // Copy the mxData array and terminate:
  memcpy(W, mxGetData(S), Len * sizeof(mxChar));
  W[Len] = L'\0';

  // mxChar to 4 byte wchar_t under Linux (wchar_t seems to be unsual here!):
  //   mxChar *M = (mxChar *) mxGetData(S);
  //   for (i = 0; i < Len; i++) {
  //      W[i] = (wchar_t) M[i];
  //   }
  
#else  // UTF8 for Linux and MacOS:
  W = MX_ARRAY_TO_STRING(S);
  if (W == NULL) {
     mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
                       ERR_HEAD "Cannot get memory for unicode string.");
  }
  
#endif
  
  return W;
}

// =============================================================================
#if WINDOWS_API_METHOD
void Core(wchar_t *FileName, int64_T Distance, Mode_t Mode, int *StatusP,
          const char **ErrMsgP)
{
  // Windows API method: SetFilePointerEx and SetEndOfFile
  
  HANDLE hFile;
  BOOL   Success;
  LARGE_INTEGER liDistance;
  DWORD  Attrib;
  
  liDistance.QuadPart = Distance;
  
  // Open the file for writing:
  hFile = CreateFileW(
             (LPCWSTR) FileName,     // pointer to name of the file
             GENERIC_WRITE,          // access (read-write) mode
             0, NULL,                // share mode and security
             OPEN_EXISTING,          // do not create
             FILE_ATTRIBUTE_NORMAL,  // attributes
             NULL);                  // attribute template handle
  
  // Handle problems from opening the file:
  if (hFile == INVALID_HANDLE_VALUE) {                // Error at opening:
     // Test existence:
     Attrib = GetFileAttributesW((LPCWSTR) FileName);
     if ((Attrib == INVALID_FILE_ATTRIBUTES)) {       // Missing file or folder:
        *ErrMsgP = MSG_no_file;
        *StatusP = ID_no_file;
     } else if (Attrib & FILE_ATTRIBUTE_DIRECTORY) {  // File is a folder:
        *ErrMsgP = MSG_bad_file;
        *StatusP = ID_bad_file;
     } else {                                         // Locked file:
        *ErrMsgP = MSG_locked;
        *StatusP = ID_locked;
     }
     return;
  }
  
  // Move the file pointer:
  if (Mode == MODE_BEGIN) {
     Success = SetFilePointerEx(hFile, liDistance, NULL, FILE_BEGIN);
  } else {
     Success = SetFilePointerEx(hFile, liDistance, NULL, FILE_END);
  }
  
  // Set the end of file:
  if (Success) {        // Movong the file pointer worked:
     Success = SetEndOfFile(hFile);
     if (!Success) {    // Setting the end-of-file failed:
        *StatusP = ID_disk_full;
        *ErrMsgP = MSG_disk_full;
     }
  } else {              // Moving the file pointer failed:
     *StatusP = ID_neg_size;
     *ErrMsgP = MSG_neg_size;
  }
  
  // Close the file:
  CloseHandle(hFile);
    
  return;
}

// =============================================================================
#elif WINDOWS_CRT_METHOD
void Core(XCHAR *FileName, int64_T Distance, Mode_t Mode, int *StatusP,
          const char **ErrMsgP)
{
  // Windows CRT method to change the file size.
  
  int     fd;
  errno_t err;
  
  fd = _wopen(FileName, _O_BINARY | _O_RDWR);
  if (fd == -1) {
     fd = _wopen(FileName, _O_BINARY | _O_RDONLY);
     if (fd == -1) {
       *ErrMsgP = MSG_locked;
       *StatusP = ID_locked;
     } else {
       *ErrMsgP = MSG_bad_file;
       *StatusP = ID_bad_file;
       _close(fd);
     }
     return;
  }
  
  if (Mode == MODE_END) {  // Move relative to the end of the file:
     Distance += _filelengthi64(fd);
  }
  
  if (Distance < 0) {      // _chsize_s crashs on negative file size
     *StatusP = ID_neg_size;
     *ErrMsgP = MSG_neg_size;
     _close(fd);
     return;
  }
  
  // Change the size and close the file (_chsize limited to 32 bit file sizes):
  err = _chsize_s(fd, Distance);
  
  _close(fd);
  
  switch (err) {   //  _chsize_s failed:
     case 0:
        break;
     case EBADF:   // Fall through, write-protection caught above already
     case EACCES:  // Actually checked above already
        *StatusP = ID_bad_file;
        *ErrMsgP = MSG_bad_file;
        break;
     case ENOSPC:
        *StatusP = ID_disk_full;
        *ErrMsgP = MSG_disk_full;
        break;
     case EINVAL:
        *StatusP = ID_neg_size;
        *ErrMsgP = MSG_neg_size;
        break;
     default:
        *StatusP = ID_unknown;
        *ErrMsgP = MSG_unknown;
  }
  
  return;
}

// =============================================================================
#elif POSIX_METHOD
void Core(XCHAR *FileName, int64_T Distance, Mode_t Mode, int *StatusP,
          const char **ErrMsgP)
{
  // Linux, MacOS or Windows CRT method to change the file size.
  
  int        fd;
  structStat statbuf;
  
  errno = 0;
  
  if (Mode == MODE_END) {  // Move relative to the end of the file:
     if (getFileStat(filename, &statbuf) == 0) {
        Distance += statbuf.st_size;
     } else {              // File or folder not existing:
        *StatusP = ID_no_file;
        *ErrMsgP = MSG_no_file;
        return;
     }
  }
  
  // Truncate or expand:
  if (truncate64(FileName, Distance) != 0) {
     // Error handling:
     switch (errno) {
        case EACCES:
           *StatusP = ID_locked;
           *ErrMsgP = MSG_locked;
           break;
        case EROFS:   // Fall through
        case EISDIR:
           *StatusP = ID_bad_file;
           *ErrMsgP = MSG_bad_file;
           break;
        case ENOENT:
           *StatusP = ID_no_file;
           *ErrMsgP = MSG_no_file;
           break;
        case EINVAL:
           *StatusP = ID_neg_size;
           *ErrMsgP = MSG_neg_size;
           break;
        default:
           *StatusP = ID_unknown;
           *ErrMsgP = MSG_unknown;
     }
  }
  
  return;
}

#else
#error No Method specified. Enable WINDOWS_API, WINDOWS_CRT or POSIX method!
#endif
