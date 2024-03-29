function red=moveRedForwardOrBack(red,red_start_to_end,move,moveOrDrop)

f=fieldnames(red);

switch moveOrDrop
    case 'move'
        MoveTheseInds=red_start_to_end;
        MoveToInds=red_start_to_end+move;
        if MoveToInds(1)<1
            MoveTheseInds(1)=abs(1-MoveToInds(1))+1;
            MoveToInds(1)=1;
        end
        if MoveToInds(2)>size(red.(f{1}),1)
            MoveTheseInds(2)=MoveTheseInds(2)-abs(size(red.(f{1}),1)-MoveToInds(2));
            MoveToInds(2)=size(red.(f{1}),1);
        end
        for i=1:size(f)
            temp=red.(f{i});
            temp(MoveToInds(1):MoveToInds(2),:)=temp(MoveTheseInds(1):MoveTheseInds(2),:);
            red.(f{i})=temp;
        end
    case 'drop'
        disp('Dropping all BUT the specified index range');
        for i=1:size(f)
            temp=red.(f{i});
            if length(temp)==1
                red.(f{i})=temp;
                continue
            end
            if length(red_start_to_end)>2
                temp=temp(red_start_to_end,:);
            else
                temp=temp(red_start_to_end(1):red_start_to_end(2),:);
            end
            red.(f{i})=temp;
        end
end