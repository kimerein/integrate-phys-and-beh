function makeFig5figs(tensor_tomatchcuedsuccess,cued_success_Response)

% all neurons, unsorted
figure();
imagesc(tensor_tomatchcuedsuccess(:,:,1));
colormap(flipud(colormap('gray')));
title('cued success all neurons unsorted');
figure();
imagesc(tensor_tomatchcuedsuccess(:,:,2));
colormap(flipud(colormap('gray')));
title('cued failure all neurons unsorted');
figure();
imagesc(tensor_tomatchcuedsuccess(:,:,3));
colormap(flipud(colormap('gray')));
title('uncued success all neurons unsorted');
figure();
imagesc(tensor_tomatchcuedsuccess(:,:,4));
colormap(flipud(colormap('gray')));
title('uncued failure all neurons unsorted');

% 

end