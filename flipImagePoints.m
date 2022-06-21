function flipped_imgpoints=flipImagePoints(imgpoints,image_xsize,image_ysize,flipwhichway)

% assume rows of imgpoints are different points, cols are x then y

flipped_imgpoints=nan(size(imgpoints));
switch flipwhichway
    case 'horizontal'
        for i=1:size(imgpoints,1)
            flipped_imgpoints(i,:)=[image_xsize-imgpoints(i,1), imgpoints(i,2)];
        end
    case 'vertical'
        for i=1:size(imgpoints,1)
            flipped_imgpoints(i,:)=[imgpoints(i,1), image_ysize-imgpoints(i,2)];
        end
end

end