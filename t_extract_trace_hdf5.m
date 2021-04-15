 fileDir='/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/615883/03122021/145Hz/615883_145Hz00001/motion_corrected'; 
 
 %enter the directory of the ROI_list.mat file and .tif files


    

        cd( fileDir);
        
       file= uigetfile('*.mat');
       load(file)
      
        %r_out=extract_trace_hdf5(ROI_list, 1,0,0,0);
        r_out=extract_trace_hdf5(ans, 1,0,0,0);
        save('615883_03122021_145Hz_trace.mat','r_out')
        %enter the directory and name of where you want to save the result
 
    
function r_out=extract_trace_hdf5(r_in, GUIpick, rmfilt, rmgreen, isLabel)
    % LATEST VERSION 5.5.20 for Rebecca
    %GUI pick is a logical (0 or 1).  1 if you'd like to select files via
    %GUI, or 0 if it should autoselect *.tif files with m_ as a prefix in 
    %the same working directory as the tifs to use.
    %
    %rmfilt is a logical (0 or 1).  Select 1 if you want to remove all .tif
    %files that have the text "f_" in their name.
    %
    %rmgreen is a logical (0 or 1).  Select 1 if you want to remove all
    %.tif files that have the name "green" in them.
    %
    %isLabel is a logical (0 or 1).  Select 1 if your input ROI structure
    %has the field name "isLabel" for labelling cells in a particular way.
    
    if GUIpick
        [selected_files,selected_dir] = uigetfile('*.tif','MultiSelect','on'); % changed to .hdf5
    else
        selected_dir = pwd;
        if rmfilt
            filt_struct = dir('*f_*.tif'); %Remove Filtered Videos from extraction
            all_struct = dir('m_*.tif');
            selected_struct = all_struct;
            for idx = 1:numel(filt_struct)
                match = strcmp({selected_struct.name}, filt_struct(idx).name)==1;
                selected_struct(match) = [];
            end
        else
            selected_struct = dir('m_*.hdf5');
        end
        if rmgreen
            green_struct = dir('*green*.tif'); %Remove green channels from extraction
            for idx = 1:numel(green_struct)
                match = strcmp({selected_struct.name}, green_struct(idx).name)==1;
                selected_struct(match) = [];
            end
        end
        temp_cell = struct2cell(selected_struct);
        selected_files = temp_cell(1,:);
    end
    
    whole_tic = tic;
    
    if class(selected_files)=='char'
        file_list(1).name = fullfile(selected_dir,selected_files);
    else
        file_list = cell2struct(fullfile(selected_dir,selected_files),'name',1);
    end
    
    %overall BG trace
    a=r_in(1).perimeter;
    [s1,s2]=size(a);
    b=[];
    I=zeros(s1,s2);
    for i=1:numel(r_in)
        b=[b;r_in(i).pixel_idx];
        I(r_in(i).pixel_idx)=1;
    end
    c=setdiff(1:s1*s2,b);
    I=~I;
    
    
    
    
    
    
    
    
    for file_idx=1:length(file_list)
        file_idx
        
        
        filename = file_list(file_idx).name;
        fprintf(['Processing ',filename,'....\n']);
        
        suffix = strsplit(filename,'.');
        suffix = suffix{end};
        
                   
        if suffix == "tif"
           
            InfoImage = imfinfo(filename);
            NumberImages=length(InfoImage);
            f_matrix = zeros(InfoImage(1).Height,InfoImage(1).Width,NumberImages,'uint16');
           
            for i=1:NumberImages
                
                f_matrix(:,:,i) = imread(filename,'Index',i,'Info',InfoImage);
            end
            
        else
            f_matrix = h5read(filename,'/image_data');
            f_matrix = permute(f_matrix, [2,1,3]); % add this line in to extract traces from SCC motion corrected videos
            InfoImage = struct('Height',size(f_matrix,1),'Width',size(f_matrix,2));
            NumberImages = size(f_matrix,3);
            
        end
        
        f_matrix = double(reshape(f_matrix,InfoImage(1).Height*InfoImage(1).Width,NumberImages));
        
        %overall BG tace
        OBG_mask=zeros(1,InfoImage(1).Height*InfoImage(1).Width);
        OBG_mask(c) = 1;
        OBG_trace = (OBG_mask*f_matrix)/sum(OBG_mask);
        
        
        
        
        
        for roi_idx=1:numel(r_in)
            roi_idx
            current_mask = zeros(1,InfoImage(1).Height*InfoImage(1).Width);
            try
                current_mask(r_in(roi_idx).pixel_idx) = 1;
                r_out(roi_idx).pixel_idx = r_in(roi_idx).pixel_idx;
            catch
                current_mask(r_in(roi_idx).pixel_idx) = 1;
                r_out(roi_idx).pixel_idx = r_in(roi_idx).pixel_idx;
            end
            image_mask = reshape(current_mask,InfoImage(1).Height,InfoImage(1).Width);
            %Find Centroid code from
            %(https://www.mathworks.com/matlabcentral/answers/322369-find-centroid-of-binary-image)
            
            [y, x] = ndgrid(1:InfoImage(1).Height, 1:InfoImage(1).Width);
            centroid = mean([x(logical(image_mask)), y(logical(image_mask))]);
            xcent = centroid(1);
            ycent = centroid(2);
            BG_radius = 50; %Number of Pixels for Local Background Radius.  Assume 1 neuron about 10 pixels in diameter
            BG_mask = (y - ycent).^2 + (x-xcent).^2 <= BG_radius.^2;
            onlyBG_mask = (BG_mask - image_mask).*I;
            BG_vect = reshape(onlyBG_mask, 1,InfoImage(1).Height*InfoImage(1).Width);
            r_out(roi_idx).BG_idx = find(BG_vect == 1);
            %Take Out Trace Values & Update r_out
            current_trace = (current_mask*f_matrix)/sum(current_mask);
            current_BG = (BG_vect*f_matrix)/sum(BG_vect);
            r_out(roi_idx).file(file_idx).filename = filename;
            r_out(roi_idx).file(file_idx).trace = current_trace;
            r_out(roi_idx).file(file_idx).BGtrace = current_BG;
          
            
            if file_idx==1
                r_out(roi_idx).trace = current_trace;
                r_out(roi_idx).BGtrace = current_BG;
            else
                r_out(roi_idx).trace = cat(2,r_out(roi_idx).trace,current_trace);
                r_out(roi_idx).BGtrace = cat(2,r_out(roi_idx).BGtrace,current_BG);
            end
            
            %Include Label
            if isLabel
                r_out(roi_idx).isLabel = r_in(roi_idx).isLabel;
            end
            
        end
        if file_idx==1
            r_out(1).overall_BG_trace=OBG_trace;
        else
             r_out(1).overall_BG_trace=cat(2,r_out(1).overall_BG_trace,OBG_trace);
        end
        
    end
    
    for roi_idx=1:numel(r_in)
        r_out(roi_idx).color = rand(1,3);
    end
        
    fprintf(['Total loading time: ',num2str(round(toc(whole_tic),2)),' seconds.\n']);
    
end