% This function presents the probabilistic associative learning (PAL) task 
% from the EMOPRED project. (c) Irene Sophia Plank, 
% irene.plank@med.uni-muenchen.de
% Participants are asked to choose between two emotion: happiness and fear.
% The images are presented for 150ms. Before the images appear, a tone is
% played for 300ms which can be predictive of the outcome of the emotion
% recognition task that follows. There are 336 trials.
% This function takes two inputs: 
% * subID        :  a string array with the subjects PID
% * eyeTracking  :  0 or 1 to choose eye tracking or not
% The function needs Psychtoolbox to function properly. If eye tracking has
% been chosen, then a LiveTracking Eye Tracker has to be connected and the
% LiveTrackToolbox for MATLAB has to be installed. The function is meant to
% be used without an external monitor. Using a second monitor can seriously
% mess up the timing. 
% The function continually saves behavioural data as well as eye tracking
% data if that option is chosen. Both files will be placed in the "Data"
% folder. 
function taskPAL

% Get all relevant inputs. 
inVar = inputdlg({'Enter PID: ', 'Eye Tracking? 0 = no, 1 = yes', 'Positive key: '}, 'Input variables', [1 45]);
subID = convertCharsToStrings(inVar{1});
eyeTracking = str2double(inVar{2});
key_pos = inVar{3};

% Get the path of this function. 
path_src = fileparts(mfilename("fullpath")); 
idx  = strfind(path_src, filesep);
path_dat = path_src(1:(idx(end)-1));

% Initialising
fx    = 40;   % size of fixation cross     
tdur  = 0.3;  % fixation cross presentation between sound and picture
pdur  = 0.15; % picture presentation

% csv file
fid = fopen(path_dat + "\Data\PAL-BV-" + subID + "_" + datestr(datetime(),'yyyymmdd-HHMM') + ".csv", 'w');
fprintf(fid, 'subID,key_pos,trl,ut,noise,block,tone,emotion,pic_num,pic_dur,key,rt\n');

% Add eye tracking stuff, if eyeTracking is set to 1. 
if eyeTracking
    % Initialise LiveTrack
    crsLiveTrackInit;
    % Open a data file to write the data to.
    crsLiveTrackSetDataFilename(char(path_dat + "\Data\PAL-ET-" + subID + "_" + datestr(datetime(),'yyyymmdd-HHMM') + ".csv"));
end

% Read in table with stimulus order. 
tbl = readtable(path_src + "\PAL_scheme-pic.csv");

% Create high and low tones.
tone{1} = [MakeBeep(440,tdur); MakeBeep(440,tdur)];
tone{2} = [MakeBeep(660,tdur); MakeBeep(660,tdur)];

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% This use of the subfunction 'UnifyKeyNames' within KbName()
% sets the names of the keys to be the same across operating systems
% (This is useful for making our experiment compatible across computers):
KbName('UnifyKeyNames');

% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer.
% For help see: Screen Screens?
screens = Screen('Screens');

% Adjust tests of the system:
Screen('Preference','SkipSyncTests', 0);
Screen('Preference','SyncTestSettings', 0.0025);

% Draw we select the maximum of these numbers. So in a situation where we
% have two screens attached to our monitor we will draw to the external
% screen. When only one screen is attached to the monitor we will draw to
% this.
% For help see: help max
screenNumber = max(screens);

% Define grey and black.
grey  = [0.43 0.43 0.43];
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);

% And the function ListenChar(), with the single number input 2,
% stops the keyboard from sending text to the Matlab window.
ListenChar(2);
% To switch the keyboard input back on later we will use ListenChar(1).

% First of all we should set up sound in general.
% This is just a single command from PsychToolBox.
% This may produce some output in your command window,
% in which PsychToolBox checks your sound card.
InitializePsychSound;

% As before we start a 'try' block, in which we watch out for errors.
try

    % Open an on screen window and color it black
    % For help see: Screen OpenWindow?
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);
    % Open a sound device
    sound_id = PsychPortAudio('Open');

    % Hide the mouse cursor
    HideCursor(window);
    
    % Get the centre coordinate of the window in pixels
    % For help see: help RectCenter
    [xCenter, yCenter] = RectCenter(windowRect);

    % Set up alpha-blending for smooth (anti-aliased) lines
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

    % Here we set the size of the arms of our fixation cross
    fixCrossDimPix = fx;
    
    % Now we set the coordinates (these are all relative to zero we will let
    % the drawing routine center the cross in the center of our monitor for us)
    xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
    yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
    allCoords = [xCoords; yCoords];
    
    % Set the line width for our fixation cross
    lineWidthPix = round(fx/10);

    % Load all pictures
    pic_files = dir([path_src '\Stimuli\SHINEd_*']);
    pics = nan(height(tbl),1);
    for i = 1:length(pic_files)
        pic = imread([pic_files(i).folder filesep pic_files(i).name]);
        pics(i) = Screen('MakeTexture', window, pic);
    end

    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 80);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Leertaste', 'center', 'center', white);
    Screen('Flip', window);

    % Start with j and send to Emotiv
    pressed = 0;
    while pressed == 0      
        % Inside this loop, we use KbCheck to check the keyboard.
        % The first output of KbCheck is just a 0 or 1
        % to say whether or not a key has been pressed.
        % The second output is the time of the key press.
        % And the third output is a vector of 0 and 1 to say which key(s)
        % was/were pressed.
        if ~pressed
            [pressed,~,keys] = KbCheck;
        end 
        
    end
    
    % Draw the fixation cross in black, set it to the center of our screen and
    % set good quality antialiasing
    Screen('DrawLines', window, allCoords,...
        lineWidthPix, black, [xCenter yCenter], 2);

    % Flip the fixation cross to the screen
    Screen('Flip', window);
    WaitSecs(tbl.iti(1));

    if eyeTracking
        % Start streaming calibrated results
        crsLiveTrackSetResultsTypeCalibrated;
        % Start tracking
        crsLiveTrackStartTracking;
    end
    
    % Go through the trials
    for j = 1:height(tbl) % trial loop

        % tone presentation
        %%should only play one sound - check table; not play them
        %%simultaneously 
        %%first fixation cross, then sound, wait, then pic
        %%maybe give tone a number too, like we did for the pictures? (like
        %%a separate column

        % This can go into the sound card with the subfunction 'FillBuffer'.
        % The id number of the sound device is the first input, then the sound wave.
        PsychPortAudio('FillBuffer', sound_id, tone{tbl.tone_num(j)});

        % Just as with Screen(), sending data to the graphics/sound card
        % does not dislpay/play it yet.
        % For this, we need the final command to play the sound.
        % This is equivalent to Screen('Flip') for the graphics card.
        PsychPortAudio('Start', sound_id);
        if eyeTracking
            % Add a comment/trigger to the eye tracking data. 
            crsLiveTrackSetDataComment(sprintf('ton_%i_%i_%i_%i_%i_%s',...
                j,tbl.ut(j),tbl.noise(j),tbl.iti(j)*1000,tbl.tone_num(j),tbl.emotion{j}));
        end
        % Wait for tdur between tone and picture.
        WaitSecs(tdur);
 
        % Draw pic stimuli to the screen. 
        Screen('DrawTexture', window, pics(tbl.pic_num(j)));
    
        % Flip to the screen.
        t_pic = Screen('Flip', window);
        if eyeTracking
            % Add a comment/trigger to the eye tracking data. 
            crsLiveTrackSetDataComment(sprintf('pic_%i_%i_%i_%i_%i_%s',...
                j,tbl.ut(j),tbl.noise(j),tbl.iti(j)*1000,tbl.tone_num(j),tbl.emotion{j}));
        end
        WaitSecs(pdur); % picture presentation
        pic_dur = (GetSecs - t_pic) * 1000;

        % Draw the fixation cross in black, set it to the center of our screen and
        % set good quality antialiasing
        Screen('DrawLines', window, allCoords,...
            lineWidthPix, black, [xCenter yCenter], 2);

        % Flip the fixation cross to the screen
        start = Screen('Flip', window);
        if eyeTracking
            % Add a comment/trigger to the eye tracking data. 
            crsLiveTrackSetDataComment(sprintf('fix_%i_%i_%i_%.2f_%i_%s',...
                j,tbl.ut(j),tbl.noise(j),tbl.iti(j),tbl.tone_num(j),tbl.emotion{j}));
        end

        % Now we start a loop for the duration of the ITI. The first key
        % that is pressed is logged. 
        dur = 0;
        pressed = 0;
        while dur < tbl.iti(j)
            % Get for how long it has been presented
            dur = GetSecs - start;
            
            % Inside this loop, we use KbCheck to check the keyboard.
            % The first output of KbCheck is just a 0 or 1
            % to say whether or not a key has been pressed.
            % The second output is the time of the key press.
            % And the third output is a vector of 0 and 1 to say which key(s)
            % was/were pressed.
            if ~pressed
                [pressed,t_press,keys] = KbCheck;
            end 
            
        end

        % The vector of keys can be input to KbName()
        % to get the name of the key that was pressed.
        key = KbName(keys);
        if iscell(key)
            key = key{end};
        end
        
        % First we check whether the key pressed was the escape key.
        if strcmp(key,'ESCAPE')
            
            % If it was, we generate an error, which will stop our script.
            % This will also send us to the catch section below,
            % where we close the screen and re-enable the keyboard for Matlab.
            error('Experiment aborted.')

        elseif strcmp(key,'NumLock')

            % If NumLock is sent, we abort because of problems with the
            % number pad.
            error('Experiment aborted: PROBLEM WITH NUMBER PAD.')
            
        end
        
        % calculate the participant's reaction time.
        rt = (t_press - t_pic) * 1000; % to get ms

        % Log all information in the log file
        fprintf(fid, '%s,%s,%i,%i,%i,%i,%s,%s,%i,%.0f,%s,%.0f\n',subID,key_pos,j,tbl.ut(j),...
            tbl.noise(j),tbl.block(j),tbl.tone{j},tbl.emotion{j},...
            tbl.pic_num(j),pic_dur,key,rt);

    end
    
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 80);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Leertaste', 'center', 'center', white);
    Screen('Flip', window);

    % End with n and send to Emotiv
    pressed = 0;
    while pressed == 0      
        % Inside this loop, we use KbCheck to check the keyboard.
        % The first output of KbCheck is just a 0 or 1
        % to say whether or not a key has been pressed.
        % The second output is the time of the key press.
        % And the third output is a vector of 0 and 1 to say which key(s)
        % was/were pressed.
        if ~pressed
            [pressed,~,~] = KbCheck;
        end 
    end
    
    % If we encounter an error...
    catch my_error

    % Close all open objects
    Screen('CloseAll');
    PsychPortAudio('Close', sound_id);

    % Show the mouse cursor
    if exist('window', 'var')
        ShowCursor(window);
    end
    
    % Clear the screen (so we can see what we are doing).
    sca;
    
    % In addition, re-enable keyboard input (so we can type).
    ListenChar(1);

    % Close the open csv file 
    fclose(fid);

    % Stop eye tracking. 
    if eyeTracking
        crsLiveTrackStopTracking;
        crsLiveTrackCloseDataFile;
        crsLiveTrackClose;
    end
    
    % Tell us what the error was.
    rethrow(my_error)
    
end

% Draw the fixation cross in black to have for two seconds before ending
Screen('DrawLines', window, allCoords,...
    lineWidthPix, black, [xCenter yCenter], 2);
Screen('Flip', window);
WaitSecs(2);

% At the end, clear the screen and re-enable the keyboard.
Screen('CloseAll');
PsychPortAudio('Close', sound_id);
ShowCursor(window);
sca;
ListenChar(1);
fclose(fid);
if eyeTracking
    crsLiveTrackStopTracking;
    crsLiveTrackCloseDataFile;
    crsLiveTrackClose;
end

end
