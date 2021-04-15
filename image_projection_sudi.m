function image_projection_sudi(fileDir)
    %Code to create a max-min projection of images across  images
     cd(fileDir); 
     video_matrix = tiffs2matrix([], 1, 0);
     final_image = max(video_matrix,[],3) - min(video_matrix,[],3);
        
     imwrite(final_image,'full_projection_m_f_max_min_batch.tif','TIFF');
     close all;
end


function [f_matrix] = tiffs2matrix(keyword, method, onlyLast5)

    if isempty(keyword)
        file_list = uigetfile('*.tif','MultiSelect','on');
        file_list = cell2struct(file_list,'name',1);
    else
        file_list = dir([keyword,'*.*']);
    end
    
    if onlyLast5 %Only Process last 5 movies to make total video
       n_files = numel(file_list);
       start_file = n_files-4;
       file_list = file_list(start_file:end);
    end

    frame_number = zeros((length(file_list)+1),1);

    InfoImage=0;
    for f=1:length(file_list)
        clear InfoImage
        filename = file_list(f).name;
        fprintf(['preparing ',filename,'\n']);
        InfoImage = imfinfo(filename);
        frame_number(f+1) = length(InfoImage);
    end

    total_frame = sum(frame_number);

    f_matrix = zeros(InfoImage(1).Height,InfoImage(1).Width,total_frame,'uint16');

    InfoImage=0;
    switch method
        case 1
            for f=1:length(file_list)
                
                clear InfoImage
                filename = file_list(f).name;
                fprintf(['Loading ',filename,'\n']);
                InfoImage = imfinfo(filename);
                NumberImages=length(InfoImage);
                frame_zero = sum(frame_number(1:f));
                for i=1:NumberImages
                    f_matrix(:,:,i+frame_zero) = imread(filename,'Index',i,'Info',InfoImage);
                end
            end
        case 2
            warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
            for f=1:length(file_list)
                clear InfoImage
                filename = file_list(f).name;
                InfoImage = imfinfo(filename);
                NumberImages=length(InfoImage);
                frame_zero = sum(frame_number(1:f));
                TifLink = Tiff(filename, 'r');
                for i = 1:NumberImages
                    TifLink.setDirectory(i);
                    f_matrix(:,:,i+frame_zero) = TifLink.read();
                end
            end
            warning('on','all');

    end
    
end

function roi_overlap(roi_list, video_matrix)

    if isempty(roi_list)
        [roi_file,roi_filedir] = uigetfile(['*.mat'],'MultiSelect','off');
        load(fullfile(roi_filedir,roi_file));
        if exist('roi_list','var')
            roi_list = r_out;
        end
    end
    
    if isempty(video_matrix)
        [video_file_list_temp,video_filedir] = uigetfile(['*.tif'],'MultiSelect','on');
        switch class(video_file_list_temp)
            case 'char'
                video_file_list{1} = video_file_list_temp;
            case 'cell'
                video_file_list = video_file_list_temp;
        end
        
        
        frame_number = zeros((length(video_file_list)+1),1);

        for video_file_idx=1:length(video_file_list)
            filename = fullfile(video_filedir,video_file_list{video_file_idx});
            InfoImage = imfinfo(filename);
            frame_number(video_file_idx+1) = length(InfoImage);
        end

        total_frame = sum(frame_number);

        video_matrix = zeros(InfoImage(1).Height,InfoImage(1).Width,total_frame,'uint16');
        
        for video_file_idx=1:length(video_file_list)
            filename = fullfile(video_filedir,video_file_list{video_file_idx});
            InfoImage = imfinfo(filename);
            NumberImages = length(InfoImage);
            frame_zero = sum(frame_number(1:video_file_idx));
            for i=1:NumberImages
                video_matrix(:,:,i+frame_zero) = imread(filename,'Index',i,'Info',InfoImage);
            end
        end
        
    end
    
    projected_image = max(video_matrix,[],3);
    
    final_image = repmat(projected_image,[1,1,3]);
    intensity_range = max(final_image(:))-min(final_image(:));
    se = strel('disk', 1, 0);
    addcolor = ~isfield(roi_list, 'color');
    
    for roi_idx=1:numel(roi_list)
        roi_mask = zeros(1024,1024);
        roi_mask(roi_list(roi_idx).PixelIdxList)=1;
        roi_boundary = imdilate(roi_mask,se)-roi_mask;
        final_image(logical(roi_boundary)) = 0;
        if addcolor
            roi_list(roi_idx).color = rand(1,3);
        end
        for color_idx=1:3
            color_image = final_image(:,:,color_idx);
            color_image(logical(roi_boundary)) = roi_list(roi_idx).color(color_idx)*intensity_range;
            final_image(:,:,color_idx) = color_image;
        end     
    end
    
    fig_1 = figure;
    imagesc(final_image);
    axis image;
    saveas(fig_1,'roi_overlay_circle.fig');
    close all;

end








