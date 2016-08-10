% Resize an existing file
%
% [Status, Msg] = FileResize(File, Bytes, [Mode])
% INPUT:
%   File:  Name of an existing file as string, with or without absolute or
%          relative path. Unicode and UNC paths are accepted.
%   Bytes: Number of bytes as numerical scalar.
%          Use UINT64, if the file size exceeds 9.01 PetaBytes.  :-)
%   Mode:  String to describe the mode:
%          'set':  Set the file size to the absolute size.
%          'move': Append or crop bytes relative to the end. If [Bytes] is
%                  negative, the file size is reduced.
%          Optional, default: 'set'.
%
% OUTPUT:
%   Status: Scalar DOUBLE. Optional, if the output is not caught by the caller,
%           this function stops with an error on problems.
%            0: Success
%           -1: Path is a folder or file is write protected
%           -2: File not found
%           -3: No space on device
%           -4: File length is negative
%           -5: File is write protected or in use
%           -6: Unknown proble
%   Msg: String, empty on success, some information in case of problems.
%
% NOTES:
% - Appended bytes are filled by zeros.
% - Resizing a file is more efficient than appending zeros manually.
%
% EXAMPLES:
%   File = tempname; fid = fopen(File, 'w'); fclose(fid);
%   FileResize(File,  8);           % 8 zeros
%   FileResize(File,  8, 'move');   % 16 zeros
%   FileResize(File, -6, 'move');   % 10 zeros
%   delete(File);
%
% COMPILE:
% This function mus be compiled before it can be used. See FileResize.c
%
% TEST:
% Run the unit-test function uTest_FileResize to check validity any speed.
%
% Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
% Author: Jan Simon, Heidelberg, (C) 2012 matlab.THISYEAR(a)nMINUSsimon.de

% $JRev: R-c V:002 Sum:ydOmQ/shKeJO Date:03-Jul-2012 02:08:49 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_FileResize $
% $File: Tools\GLFile\FileResize.m $
% History:
% 001: 03-Jul-2012 00:58, First version.
