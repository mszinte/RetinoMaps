function [expDes] = runTrials(scr,const,expDes,my_key)
% ----------------------------------------------------------------------
% [expDes]=runTrials(scr,const,expDes,my_key)
% ----------------------------------------------------------------------
% Goal of the function :
% Draw stimuli of each indivual trial and waiting for inputs
% ----------------------------------------------------------------------
% Input(s) :
% scr : struct containing screen configurations
% const : struct containing constant configurations
% expDes : struct containg experimental design
% my_key : structure containing keyboard configurations
% ----------------------------------------------------------------------
% Output(s):
% resMat : experimental results (see below)
% expDes : struct containing all the variable design configurations.
% ----------------------------------------------------------------------
% Function created by Martin SZINTE (martin.szinte@gmail.com)
% Last update : 22 / 09 / 2020
% Project :     pRFexp
% Version :     1.0
% ----------------------------------------------------------------------

for bar_pass = 1:const.bar_dir_num
    
    % Write in log/edf
    log_txt     =   sprintf('bar pass %i started',bar_pass);
    
    if const.tracker
        Eyelink('message','%s',log_txt);
    end
    
    % Compute and simplify var and rand
    % ---------------------------------

    % trials number in this bar pass    
    bar_trials              =   expDes.expMat(:,5) == bar_pass;
    bar_trials_num          =   expDes.expMat(bar_trials,2);

    % Cond1 : Task
    cond1                   =   expDes.expMat(bar_trials,3);

    % Var 1 : Bar direction
    var1                    =   expDes.expMat(bar_trials,4);

    % Rand 1 : Stimulus orientation
    rand1                   =   expDes.expMat(bar_trials,7);

    if const.checkTrial && const.expStart == 0
        fprintf(1,'\n\n\bar_pass========================  BAR PASS %3.0f ========================\n',bar_pass);
        fprintf(1,'\n\tTask                         =\bar_pass');
        fprintf(1,'%s\bar_pass',expDes.txt_cond1{cond1(1)});
        fprintf(1,'\n\tBar direction                =\bar_pass');
        fprintf(1,'%s\bar_pass',expDes.txt_var1{var1(1)});
        fprintf(1,'\n\tStimulus orientation         =\bar_pass');
        fprintf(1,'%s\bar_pass',expDes.txt_rand1{rand1});
    end

    % Prepare stimuli
    % ---------------
    % Stimulus

    if var1(1) == 1 || var1(1) == 5 || var1(1) == 9
        rand_num_tex_cond       =   const.rand_num_tex_hor;
        bar_step_val_cond       =   1:const.bar_step_hor;
        pre_txt_cond            =   'hor_';
        bar_step_cond           =   const.bar_steps_hor;
        time2probe_cond         =   const.time2probe_hor;
        time2resp_cond          =   const.time2resp_hor;
        resp_reset_cond         =   const.resp_reset_hor;
        trial_start_cond        =   const.trial_start_hor;
        trial_end_cond          =   const.trial_end_hor;
        probe_onset_cond        =   const.probe_onset_hor;

    elseif var1(1) == 3 || var1(1) == 7
        rand_num_tex_cond       =   const.rand_num_tex_ver;
        bar_step_val_cond       =   1:const.bar_step_ver;
        pre_txt_cond            =   'ver_';
        bar_step_cond           =   const.bar_steps_ver;
        time2probe_cond         =   const.time2probe_ver;
        time2resp_cond          =   const.time2resp_ver;
        resp_reset_cond         =   const.resp_reset_ver;
        trial_start_cond        =   const.trial_start_ver;
        trial_end_cond          =   const.trial_end_ver;
        probe_onset_cond        =   const.probe_onset_ver;
        
    end

    % define image of next frame
    bar_step                =   1;

    % define order of the bar (only one direction is saved as img) + duration of the trial
    bar_step_val_rev_cond   =   bar_step_val_cond(end:-1:1);
    rand_num_tex            =   rand_num_tex_cond;
    if var1(bar_step) == 1 || var1(bar_step) == 3
        drawf_max = const.num_drawf_max_hor;
    elseif var1(bar_step) == 5 || var1(bar_step) == 7
        drawf_max = const.num_drawf_max_ver;
    elseif var1(bar_step) == 9
        drawf_max = const.num_drawf_max_blk;
    end
    
    % wait for bar_pass press in trial beginning
    if bar_pass == 1
        Screen('FillRect',scr.main,const.background_color);
        drawEmptyTarget(scr,const,scr.x_mid,scr.y_mid);
        Screen('Flip',scr.main);
        first_trigger           =   0;
        expDes.mri_band_val     =   my_key.first_val(3);

        while ~first_trigger
            if const.scanner == 0 || const.scannerTest
                first_trigger           =   1;
                mri_band_val            =   -8;
                
            else
                keyPressed              =   0;
                keyCode                 =   zeros(1,my_key.keyCodeNum);
                for keyb = 1:size(my_key.keyboard_idx,2)
                    [keyP, keyC]            =   KbQueueCheck(my_key.keyboard_idx(keyb));
                    keyPressed              =   keyPressed+keyP;
                    keyCode                 =   keyCode+keyC;
                end
                if const.scanner == 1
                    input_return = [my_key.ni_session2.inputSingleScan,my_key.ni_session1.inputSingleScan];

                    if input_return(my_key.idx_mri_bands) == ~expDes.mri_band_val
                        keyPressed              = 1;
                        keyCode(my_key.mri_tr)  = 1;
                        expDes.mri_band_val     = ~expDes.mri_band_val;
                        mri_band_val            = input_return(my_key.idx_mri_bands);
                    end
                end

                if keyPressed
                    if keyCode(my_key.escape) && const.expStart == 0
                        overDone(const,my_key)
                    elseif keyCode(my_key.mri_tr)
                        first_trigger          =   1;
                    end
                end
            end
        end
        
        % write in log/edf
        log_txt                 =   sprintf('bar pass %i event mri_trigger val = %i',bar_pass,mri_band_val);
        if const.tracker
            Eyelink('message','%s',log_txt);
        end
    
        t_start = GetSecs;
        vbl = GetSecs;
        log_txt  =   sprintf('bar pass %i trial onset %i',bar_pass,bar_trials_num(bar_step));
        if const.tracker
            Eyelink('message','%s',log_txt);
        end
        expDes.expMat(bar_trials_num(bar_step),8) = vbl;
    end
    
    drawf = 1;
    while drawf <= drawf_max
         
        % get bar step value
        bar_step =   bar_step_cond(drawf,bar_pass);
        
        % load stim
        % define order of the bar (only one direction is saved as img)
        if var1(bar_step) == 1 || var1(bar_step) == 3
            bar_step_num            =   bar_step_val_cond(bar_step);
        elseif var1(bar_step) == 5 || var1(bar_step) == 7
            bar_step_num            =   bar_step_val_rev_cond(bar_step);
        end

        % define displayed direction of the bar (only 0 deg direcion motion
        % reversing the bar and fix stim reversed)
        angle                   =   0;
        stim_ori                =   rand1(bar_step);
        
        % define when to reset resp
        if resp_reset_cond(drawf,bar_pass)
            resp                =   0;
        end

        % save stim staircase level
        if probe_onset_cond(drawf,bar_pass)
            expDes.expMat(bar_trials_num(bar_step),10)  =   expDes.stim_stair_val;
        end
        
        % define name of next frame
        if var1(bar_step) == 9
            screen_filename         =   sprintf('%s/blank.mat',const.stim_folder);
        else
            if time2probe_cond(drawf,bar_pass)
                screen_filename         =   sprintf('%s/%sprobe_barStep-%i_kappaStim%i_stimOri-%i_noiseRand%i.mat',...
                    const.stim_folder,pre_txt_cond,bar_step_num,expDes.stim_stair_val,stim_ori,rand_num_tex(drawf));
            else
                screen_filename         =   sprintf('%s/%snoprobe_barStep-%i_noiseRand%i.mat',...
                    const.stim_folder,pre_txt_cond,bar_step_num,rand_num_tex(drawf));
            end
        end
        
        % load the matrix
        load(screen_filename,'screen_stim');
        
        % make texture
        expDes.tex = Screen('MakeTexture',scr.main,screen_stim,[],[],[],angle);
        
        % draw texture
        Screen('DrawTexture',scr.main,expDes.tex,[],const.stim_rect)        
        
        
        % Check keyboard
        % --------------
        keyPressed              =   0;
        keyCode                 =   zeros(1,my_key.keyCodeNum);
        for keyb = 1:size(my_key.keyboard_idx,2)
            [keyP, keyC]            =   KbQueueCheck(my_key.keyboard_idx(keyb));
            keyPressed              =   keyPressed+keyP;
            keyCode                 =   keyCode+keyC;
        end

        if const.scanner == 1
            input_return =  [my_key.ni_session2.inputSingleScan,my_key.ni_session1.inputSingleScan];

            % button press
            if input_return(my_key.idx_button_right1) == my_key.button_press_val
                keyPressed              = 1;
                keyCode(my_key.right1)  = 1;
            elseif input_return(my_key.idx_button_left4) == my_key.button_press_val
                keyPressed              =  1;
                keyCode(my_key.left4)   =  1;
            end

            % mri trigger
            if input_return(my_key.idx_mri_bands) == ~expDes.mri_band_val
                keyPressed              = 1;
                keyCode(my_key.mri_tr)  = 1;
                expDes.mri_band_val     = ~expDes.mri_band_val;
                mri_band_val            = input_return(my_key.idx_mri_bands);
            end
        end

        if keyPressed
            if keyCode(my_key.mri_tr)
                % write in log/edf
                log_txt                 =   sprintf('bar pass %i event mri_trigger val = %i',bar_pass,mri_band_val);
                if const.tracker
                    Eyelink('message','%s',log_txt);
                end
            end
            if keyCode(my_key.escape)
                if const.expStart == 0
                    overDone(const,my_key)
                end
            elseif keyCode(my_key.left4)
                % update staircase
                if time2resp_cond(drawf,bar_pass) && resp == 0
                    % write in log/edf
                    log_txt                 =   sprintf('bar pass %i trial %i event %s',bar_pass, bar_trials_num(bar_step), my_key.left4Val);
                    if const.tracker
                        Eyelink('message','%s',log_txt);
                    end

                    if cond1(bar_step) == 1
                        switch rand1(bar_step)
                            case 1; response                =   0;
                            case 2; response                =   1;
                        end
                    end
                    expDes.expMat(bar_trials_num(bar_step),11)  =   response;
                    expDes.expMat(bar_trials_num(bar_step),13)  =   GetSecs;
                    [expDes]                =   updateStaircase(cond1(bar_step),const,expDes,response);
                    resp                    =   1;
                end
            elseif keyCode(my_key.right1)  % cw button
                % update staircase
                if time2resp_cond(drawf,bar_pass) && resp == 0
                    % write in log/edf
                    log_txt                 =   sprintf('bar pass %i trial %i event %s',bar_pass,bar_trials_num(bar_step), my_key.right1Val);
                    if const.tracker
                        Eyelink('message','%s',log_txt);
                    end
                    if cond1(bar_step) == 1
                        switch rand1(bar_step)
                            case 1; response                =   1;
                            case 2; response                =   0;
                        end
                    end
                    expDes.expMat(bar_trials_num(bar_step),11)  =   response;
                    expDes.expMat(bar_trials_num(bar_step),13)  =   GetSecs;
                    [expDes]                =   updateStaircase(cond1(bar_step),const,expDes,response);
                    resp                    =   1;
                end
            end
        end
        
        % Create movie
        if const.mkVideo
            expDes.vid_num          =   expDes.vid_num + 1;
            image_vid               =   Screen('GetImage', scr.main);
            imwrite(image_vid,sprintf('%s_frame_%i.png',const.movie_image_file,expDes.vid_num))
            writeVideo(const.vid_obj,image_vid);
        end
        
        % flip screen
        when2flip = vbl + const.patch_dur - scr.frame_duration/2;
        vbl = Screen('Flip',scr.main, when2flip);

        % probe
        if probe_onset_cond(drawf,bar_pass)
            
            % probe onset
            if cond1(bar_step) == 1
                log_txt                 =   sprintf('bar pass %i stimulus probe onset %i',bar_pass,bar_trials_num(bar_step));
            end
            if const.tracker
                Eyelink('message','%s',log_txt);
            end
            expDes.expMat(bar_trials_num(bar_step),12) = vbl;
        end
        
        if trial_end_cond(drawf,bar_pass)
            % trial offset
            log_txt = sprintf('bar pass %i trial offset %i',bar_pass, bar_trials_num(bar_step));
            if const.tracker
                Eyelink('message','%s',log_txt);
            end
            expDes.expMat(bar_trials_num(bar_step),9) = vbl;
        end
        
        % trial onset
        if trial_start_cond(drawf,bar_pass) && ~(bar_pass == const.bar_dir_num && drawf == drawf_max)
            if drawf == drawf_max
                log_txt  =   sprintf('bar pass %i trial onset %i',bar_pass+1, bar_trials_num(bar_step)+1);
            else
                log_txt  =   sprintf('bar pass %i trial onset %i',bar_pass, bar_trials_num(bar_step)+1);
            end
            if const.tracker
                Eyelink('message','%s',log_txt);
            end
            
            expDes.expMat(bar_trials_num(bar_step)+1,8) = vbl;
        end
        
        % Drawing frame number
        drawf = drawf + 1;

    end
    if const.checkTrial
        fprintf(1,'num draw = %i, dur = %1.2f',drawf, GetSecs - t_start);
    end

    % write in log/edf
    log_txt                     =   sprintf('bar pass %i stopped',bar_pass);
    if const.tracker
        Eyelink('message', '%s',log_txt);
    end

end
end