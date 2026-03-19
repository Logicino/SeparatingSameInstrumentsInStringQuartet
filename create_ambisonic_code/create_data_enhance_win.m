close all
clear all
clc

%create evaluation/validation set with stereo reverb added
%using URMP and Bach10 as anechoic isolated instruments
%And openair for b-format impulse responses
%Trevor Cox, University of Salford for Cadenza project July 2024

fs = 44100;     %target sample rate

rng(0);   %set seed so each run identical

rootdir = 'F:\Cadenza\data_enhance';
dir(rootdir)
irlist = dir(fullfile(rootdir, '**\*.wav'))
filelist = dir(fullfile(rootdir, '**\*.flac'))
% filelist = dir(fullfile(rootdir, '**\*.wav'));  %get list of files and folders in any subfolder
% print(filelist)
irlist = irlist(~[irlist.isdir]); 
filelist = filelist(~[filelist.isdir]);         %remove folders from list

%where the code says "URMP" in a variable name it also works for "Bach10"
% 在所有的wav文件里分别筛选出URMP和IR的文件
Ensemblefilelist = [];                              %sort lists into solo instrument audio and impulse responses
IRfilelist = [];

% 删除以'.'开头的文件
% 1. 删除filelist
for n=1:length(filelist)
    % Check if the file name starts with a dot
    if strncmp(filelist(n).name, '.', 1)
        continue; % Skip this file
    end

    Ensemblefilelist = [Ensemblefilelist ; filelist(n)]; %URMP instruments

end

% 2. 删除irlist
for n=1:length(irlist)
    % Check if the file name starts with a dot
    if strncmp(irlist(n).name, '.', 1)
        continue; % Skip this file
    end

    IRfilelist = [IRfilelist ; irlist(n)]; %URMP instruments

end

% n = 1;                                          %there are very short files in the URMP list that are removed here from list
% while n<length(URMPfilelist)
%     if URMPfilelist(n).bytes<1000
%         URMPfilelist(n:end-1) = URMPfilelist(n+1:end);
%         URMPfilelist = URMPfilelist(1:end-1);
%     else
%         n = n + 1;
%     end
% end

% 移除mix的文件
n = 1;                                          %don't need to process mixes in the folders remove them from list
while n<length(Ensemblefilelist) + 1
    if strfind(lower(Ensemblefilelist(n).name),'mix')
        Ensemblefilelist(n:end-1) = Ensemblefilelist(n+1:end);
        Ensemblefilelist = Ensemblefilelist(1:end-1);
    else
        n = n + 1;
    end
end

%%% 至此，已经剔除所有mix和.开头文件

root = 'Stereo_Reverb_Ensembleset';                   %create directory for eval/val set

n = 1;                                  %loop over 10all the audio files
while n <= length(Ensemblefilelist)

    m = n+1;             %find out how many instruments in the ensemble by looking at number of .wav in the folder
    while m<=length(Ensemblefilelist) && ~isempty(strfind(Ensemblefilelist(n).folder,Ensemblefilelist(m).folder))
        m = m + 1;
    end
    Ninstruments = m-n;     %no of instruments in the ensemble
    % load audio
    for m=1:Ninstruments
        Ensemblefilelist(m+n-1).name
        [x1,fs1] = audioread([Ensemblefilelist(m+n-1).folder '\' Ensemblefilelist(m+n-1).name]);   % 加载音频
        if abs(fs1-fs)>1
            x1 = audioresample(x1,InputRate=fs1,OutputRate=fs);      %get everything into 44.1kHz if not already 重采样  111
        end
        x11 = x1(:,1);
        x12 = x1(:,2);
        % x = x1;
        x_mono = (x11 + x12) / 2;
        x(:,m) = x_mono;
    end
    rmsx = rms(x(:));           %calculate rms across all instruments for later normalisation

    rootout = Ensemblefilelist(n).folder;       %create output folder for the reverberated audio
    rootout = replace(rootout,'EnsembleSet',[root '\EnsembleSet']);
    % rootout = replace(rootout,'BACH10',[root '/BACH10']);
    mkdir(rootout);    %output directory

    %randommly choose impulse response 随机选取一个B-FormatIR文件
    NIR = floor(rand(1,1)*length(IRfilelist))+1;
    % load b-format impulse response ( 4 channel ) 读取IR文件
    [y,fs2] = audioread([IRfilelist(NIR).folder '\' IRfilelist(NIR).name]);
    if abs(fs2-fs)>1
        y = audioresample(y,InputRate=fs2,OutputRate=fs);      %get everything into 44.1kHz
    end
    fid = fopen([rootout '\readme.txt'],'wt');
    fprintf(fid,['Reverberation added from: ' IRfilelist(NIR).name]);
    fclose(fid);

    %apply impulse responses using Fourier
    N = max([length(x) length(y)]);     %length for fft
    mx = zeros(N,2);        %mix        %to hold mix

    %angular positions of sources
    source_angles = [1:Ninstruments]'*10;                   %angles for sources, ten degree spacing
    source_angles = source_angles - mean(source_angles);    %centre on zero

    for m=1:Ninstruments

        % bformat rotate. This is how I am positioning the instruments on
        % the stage. There is no correct way of doing this. Rotating makes
        % the direct sound in the correct direction but the room
        % reflections are wrong.
        y_rot = rotateBformat(y, source_angles(m), 0, 0);
        % convert b-format to stereo
        %s = convertBFormToStereoUHJ(y_rot);    %function from York but has
        %a n error in it
        s = convertBFormToMSStereo(y_rot);  %b-format to mid-side stereo
        S = fft(s,N);                       %get spectra of stereo impulse responses

        X = fft(x(:,m));                    %spectrum of mono instrument audio
        for nchannel = 1:2
            Z = X.*S(:,nchannel);           %convolution (in freq domain multiplication)
            z(:,nchannel,m) = ifft(Z);      %back to time domain
        end
        mx = mx + z(:,:,m);                 %make mix through accumulation

    end

    %RMS normalisation not being used
    % rmsz = rms(z(:));                            %rms across all instruments
    % z = z*rmsx/rmsz/Ninstruments;                        %normalise back to original rms /Ninstruments so mix doesn't clip
    % mx = mx*rmsx/rmsz/Ninstruments;                      %normalise back to original rms /Ninstruments so mix doesn't clip

    %Bach10 and URMP have different normalisation strategies
    %Peak(abs(mix))=1 normalisation. With instruments audio filesset so they sum to
    %make the mix.
    mxmax = max(max(abs(mx)));
    z = z/mxmax;                        %normalise so max(mix)=1
    mx = mx/mxmax;                      %normalise so max(mix)=1
    
    
    for m=1:Ninstruments                    %write audio files

        fname = [rootout '\' URMPfilelist(m+n-1).name];    %file name
        fname = replace(fname,'wav','flac');        %save output as flac. 16 bit is default

        if max(max(max(abs(z))))>1               %check for clipping issues
            display('clipping z')
            display(fname)
            pause
        end

        audiowrite(fname,z(:,:,m),fs)              %save audio

        I = strfind(URMPfilelist(n).folder,'\');    %make audio output name for mix
        fname = URMPfilelist(n).folder;
        fname = fname(I(end)+1:end);
        fname = [rootout '\mix_' fname '.flac'];

    end

    if max(max(abs(mx)))>1              %check for clipping issues
        display('clipping mx')
        display(fname)
        pause
    end

    audiowrite(fname,mx,fs)             %save audio

    n = n + Ninstruments;                   %increment counter to next set of instruments
    clear x                           %to prevent errors when these are reused 
    clear z                           %to prevent errors when these are reused 

end





