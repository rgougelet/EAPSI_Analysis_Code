% datafiles = subdir('*_data.txt');

%% Extract xml data names
dell_xmlfiles = subdir('*dell*.xml');
lt_xmlfiles = subdir('*xplane_udp.xml');

%% Test contents of xml files
for dellxml_ind = 1:length(dell_xmlfiles)
    % identify dell xml file
    [PATHSTR,NAME,EXT] = fileparts(dell_xmlfiles(dellxml_ind).name);
    % find the lt xml file in the same directory
    lt_xmlfilename = [PATHSTR,filesep,'xplane_udp.xml'];
    
    % get lt xml contents
    lt_tree = xml2struct(lt_xmlfilename);
    lt_nVars = length(lt_tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item);
    lt_varNames = cell(1,lt_nVars);
    for itemN = 1:lt_nVars
        lt_varNames{itemN} = lt_tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item{itemN}.Attributes.ref;
    end

    % get dell xml contents
    dell_tree = xml2struct(dell_xmlfiles(dellxml_ind).name);
    dell_nVars = length(dell_tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item);
    dell_varNames = cell(1,dell_nVars);
    for itemN = 1:dell_nVars
        dell_varNames{itemN} = dell_tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item{itemN}.Attributes.ref;
    end

    % compare dell vs lt xml contents, number of variables
    if dell_nVars == lt_nVars
        display('Number of variables equal.')
    else
        display('Number of variables NOT equal!!!!')
        continue
    end

    % compare dell vs lt xml contents, names of variables
    display('Comparing names of variables.')
    for ind = 1:length(lt_varNames)
        if ~strcmp(lt_varNames(ind), dell_varNames(ind))
            display([lt_xmlfiles(ltxml_ind).name, ' NOT equal'])
            continue
        end
    end
    display('All match!')
end
%go through lt xml files and compare names
%are they the exact same name?
%is the dell xml filename the same when dell_ is removed?
%compare contents of files
% indices = [];
% for ltxml_ind = 1:length(lt_xmlfiles)
%     display(['Looking at ',lt_xmlfiles(ltxml_ind).name])
%     for dellxml_ind = 1:length(dell_xmlfiles)
%         display(['    Looking at ',dell_xmlfiles(dellxml_ind).name])
%         if strcmp(lt_xmlfiles(ltxml_ind).name,dell_xmlfiles(dellxml_ind).name)
%             indices = [indices ltxml_ind];
%             continue
%         end
%         lt_tree = xml2struct(lt_xmlfiles(ltxml_ind).name); % need to download this from the web
%         lt_nVars = length(lt_tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item);
%         lt_varNames = cell(1,lt_nVars);
%         for itemN = 1:lt_nVars
%         	lt_varNames{itemN} = lt_tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item{itemN}.Attributes.ref;
%         end
%         
%         dell_tree = xml2struct(dell_xmlfiles(dellxml_ind).name); % need to download this from the web
%         dell_nVars = length(dell_tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item);
%         dell_varNames = cell(1,dell_nVars);
%         for itemN = 1:dell_nVars
%         	dell_varNames{itemN} = dell_tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item{itemN}.Attributes.ref;
%         end
%         
%         display('    Comparing number of variables.')
% %         length(lt_varNames) ~= length(dell_varNames)
%         if length(lt_varNames) ~= length(dell_varNames)
%             display([lt_xmlfiles(ltxml_ind).name, ' not equal'])
%             continue
%         end
%         display('    Equal.')
%         
%         display('    Comparing names of variables.')
%         for ind = 1:length(lt_varNames)
% %             ~strcmp(lt_varNames(ind), dell_varNames(ind))
%             if ~strcmp(lt_varNames(ind), dell_varNames(ind))
%                 display([lt_xmlfiles(ltxml_ind).name, ' not equal'])
%                 break
%             end
%         end
%         display('     All match!')
%     end
% end
% 
% % cd C:\Users\Rob\Desktop\Dropbox\EAPSI\EAPSI_MRI\MRI_Subject21\2015-07-29_15-02-26 % point to date_time folder of interest
% % tree = xml2struct('xplane_udp.xml'); % need to download this from the web
% % nVars = length(tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item);
% % varNames = cell(1,nVars);
% % for itemN = 1:nVars
% %    varNames{itemN} = tree.XdevLCoreProperties.XPlane.XPlaneDataRef.Item{itemN}.Attributes.ref;
% % end
% % 
% % %% Import data file
% % % If you get an error, go into text file and fix line with error. Seems to
% % % be problem with newline character being put in at certain times. I went 
% % % in and separated the two lines, then took the average of the last column 
% % % value of the before and after lines because it was missing. Also, delete 
% % % the last line because it will not have all of the columns.
% % dataFileName = ls('*_data.txt');
% % data = load(dataFileName);
% % 
% % %% Map variable names
% % time = data(:,1);
% % for varN = 1:nVars
% %     varName = strrep(varNames{varN}, '/', '_');
% %     eval([varName,' = data(:,1+varN);'])
% % end
% % performance = data(:,end);
% % startIndex = find(diff(sim_cockpit_radios_nav1_freq_hz(1:100)))+1;
% % startTime = time(startIndex);
% % 
% % clear data tree varN varName nItems itemN dataFileName