function uTest_FileResize(doSpeed)
% Automatic test: FileResize
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% uTest_FileResize(doSpeed)
% INPUT:
%   doSpeed: Optional logical flag to trigger time consuming speed tests.
%            Default: TRUE. If no speed tested are defined, this is ignored.
% OUTPUT:
%   On failure the test stops with an error.
%
% Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
% Author: Jan Simon, Heidelberg, (C) 2012 matlab.THISYEAR(a)nMINUSsimon.de

% $JRev: R-d V:003 Sum:IkmpUJrQwCJp Date:03-Jul-2012 01:02:35 $
% $License: BSD $
% $File: Tools\UnitTests_\uTest_FileResize.m $
% History:
% 001: 01-Jul-2012 12:59, First version.

% Initialize: ==================================================================
% Global Interface: ------------------------------------------------------------
ErrID = ['JSimon:', mfilename, ':Fault'];

% Initial values: --------------------------------------------------------------
% Program Interface: -----------------------------------------------------------
if nargin == 0
   doSpeed = true;
end

% User Interface: --------------------------------------------------------------
% Do the work: =================================================================
disp(['== Test FileResize  ', datestr(now, 0)]);
fprintf('  Version: %s\n', which('FileResize'));

File = fullfile(tempdir, 'FileResize.test');

% Cleanup file, which might exist from a former crashed test:
if FileExist(File)
   delete(File);
end

% ------------------------------------------------------------------------------
disp('== Known answer tests:');
FID = fopen(File, 'w');
if FID < 0
   error(ErrID, 'Cannot create empty file.');
end
fclose(FID);

D = dir(File);
if D.bytes ~= 0
   error(ErrID, 'Created file is not empty.');
end

try
   FileResize(File, 0);
catch
   error(ErrID, 'FileResize(0) of empty file crashed: %s', lasterr);
end

try
   FileResize(File, 0, 'move');
catch
   error(ErrID, 'FileResize(0, move) of empty file crashed: %s', lasterr);
end

try
   FileResize(File, 1);
catch
   error(ErrID, 'FileResize(1) of empty file crashed: %s', lasterr);
end
D = dir(File);
if D.bytes ~= 1
   error(ErrID, 'Size~=1 after FileResize(1) of empty file.');
end

try
   FileResize(File, 1, 'move');
catch
   error(ErrID, 'FileResize(1, move) of empty file crashed: %s', lasterr);
end
D = dir(File);
if D.bytes ~= 2
   error(ErrID, 'Size~=2 after FileResize(1, move) of 1-byte file.');
end

try
   FileResize(File, 1);
catch
   error(ErrID, 'FileResize(1) of 2-byte file crashed: %s', lasterr);
end
D = dir(File);
if D.bytes ~= 1
   error(ErrID, 'Size~=1 after FileResize(1) of 2 byte file.');
end

try
   FileResize(File, 32, 'set');
catch
   error(ErrID, 'FileResize(32, set) of 1-byte file crashed: %s', lasterr);
end
D = dir(File);
if D.bytes ~= 32
   error(ErrID, 'Size~=32 after FileResize(32, set) of 1 byte file.');
end

try
   FileResize(File, 32, 'move');
catch
   error(ErrID, 'FileResize(32, move) of 32-byte file crashed: %s', lasterr);
end
D = dir(File);
if D.bytes ~= 64
   error(ErrID, 'Size~=64 after FileResize(32,move) of 32 byte file.');
end

try
   FileResize(File, -64, 'move');
catch
   error(ErrID, 'FileResize(-64, move) of 64-byte file crashed: %s', lasterr);
end
D = dir(File);
if D.bytes ~= 0
   error(ErrID, 'Size~=0 after FileResize(-64,move) of 64 byte file.');
end

fprintf('  ok\n\n');

% Provoke problems: ------------------------------------------------------------
disp('== Provoke errors:');

% Bad trials with caught output:
try
   [Status, Msg] = FileResize(File, -1);
catch
   error(ErrID, '[S,M]=FileResize(-1) of empty file crashed: %s', lasterr);
end
if Status >= 0
   error(ErrID, '[S,M]=FileResize(-1) does not fail.');
end

tooLazy = false;
try
   FileResize(File, -1);
   tooLazy = true;
catch  % Fine
end
if tooLazy
   error(ErrID, 'FileResize(-1) did not throw an error.');
end

try
   [Status, Msg] = FileResize(File, -1, 'move'); %#ok<*NASGU>
catch
   error(ErrID, '[S,M]=FileResize(-1, move) crashed: %s', lasterr);
end
if Status >= 0
   error(ErrID, '[S,M]=FileResize(-1, move) does not fail.');
end

try
   FileResize(File, -1, 'move');
   tooLazy = true;
catch  % Fine
end
if tooLazy
   error(ErrID, 'FileResize(-1, move) did not throw an error.');
end

try
   [Status, Msg] = FileResize(File);
   tooLazy = true;
catch
end
if tooLazy
   error(ErrID, '[S,M]=FileResize(File) accepted a single input?');
end

try
   FileResize(File);
   tooLazy = true;
catch
end
if tooLazy
   error(ErrID, 'FileResize(File) accepted a single input?');
end

try
   FileResize(File, 1, 'set', 'Junk');
   tooLazy = true;
catch
end
if tooLazy
   error(ErrID, 'FileResize(File, 1, set, Junk) accepted 4 inputs?');
end

try
   [Status, Msg] = FileResize(File, 1, 'set', 'Junk');
   tooLazy = true;
catch
end
if tooLazy
   error(ErrID, '[S,M]=FileResize(File, 1, set, Junk) accepted 4 inputs?');
end

try
   [Status, Msg] = FileResize(File, 1, 'unknown');
catch
end
if Status >= 0
   error(ErrID, '[S,M]=FileResize(File, 1, unknown) accepted bad string?');
end

fileattrib(File, '-w');
try
   FileResize(File, 8);
   tooLazy = true;
catch
end
fileattrib(File, '+w');
if tooLazy
   error(ErrID, '[S,M]=FileResize(File, 8) expanded write protected file?');
end

missingFile = tempname;
try
   FileResize(missingFile, 1);
   tooLazy = true;
catch
end
if tooLazy
   error(ErrID, 'FileResize(tempname, 1) worked on missing file?');
end

% Error messages:
[Status, Msg] = FileResize(missingFile, 1);
if Status < 0
   fprintf('  Missing file detected:\n    %s\n', Msg);
else
   error(ErrID, '[S,M]=FileResize(tempname, 1) worked on missing file?');
end

[Status, Msg] = FileResize(tempdir, 1);
if Status < 0
   fprintf('  Folder detected:\n    %s\n', Msg);
else
   error(ErrID, '[S,M]=FileResize(tempdir, 1) worked on folder?');
end

[Status, Msg] = FileResize(File, -1);
if Status < 0
   fprintf('  Negative file length detected:\n    %s\n', Msg);
else
   error(ErrID, '[S,M]=FileResize(File, -1) worked with negative file length?');
end

[Status, Msg] = FileResize(File, 8);
[Status, Msg] = FileResize(File, -9, 'move');
if Status < 0
   fprintf('  Negative file length after MOVE detected:\n    %s\n', Msg);
else
   error(ErrID, ...
      '[S,M]=FileResize(File, -9, move) worked with negative file length?');
end

fileattrib(File, '-w');
[Status, Msg] = FileResize(File, 8);
fileattrib(File, '+w');
if Status < 0
   fprintf('  Write protected file detected:\n    %s\n', Msg);
else
   error(ErrID, ...
      '[S,M]=FileResize(File, 8) worked inspite of write protection?');
end

% I hesitate to write a command, which occupies more than the available disk
% space...

fprintf('  ok\n');

% Speed: -----------------------------------------------------------------------
if doSpeed
   fprintf('\n== Speed tests:\n');
   loop = 1000;
   FileResize(File, 0);
   
   for len = [10, 100, 1000, 10000];
      fprintf('  Append %5d bytes, %d loops:\n', len, loop);
      data = repmat(uint8(0), 1, len);
      
      tic;
      for i = 1:loop
         fid = fopen(File, 'a');
         fwrite(fid, data, 'uint8');
         fclose(fid);
      end
      mtime = toc;
      fprintf('    FWRITE:     %7.4f sec\n', mtime);
      
      FileResize(File, 0);
      tic;
      for i = 1:loop
         FileResize(File, len, 'move');
      end
      xtime = toc;
      
      fprintf('    FileResize: %7.4f sec\n', xtime);
   end
end

% Cleanup: ---------------------------------------------------------------------
if FileExist(File)
   delete(File);
end

fprintf('\n== FileResize has passed the tests.\n');

% ******************************************************************************
function  E = FileExist(File)
E = any(exist(File, 'file') == 2);
% return;

