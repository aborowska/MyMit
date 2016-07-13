function str = num2bank(num)
    % By  Ashkan Ziabakhshdeylami, nl.mathworks.com/matlabcentral/newsreader/view_thread/169226
    % Modyfied by Agnieszka Borowska (so that integers are still represented without decimal points)
    % str = arrayfun(@(x) num2bankScalar(x) , num, 'UniformOutput', false);
% end

% function str = num2bankScalar(num)
    num = floor(num*100)/100;
    str = num2str(num);
    k = find(str == '.', 1);

    if(isempty(k))
        str = [str,'.00'];
    end
   
        
    FIN = min(length(str),find(str == '.')-1);
    for i = FIN-2:-3:2
        str(i+1:end+1) = str(i:end);
        str(i) = ',';
    end
    
%     if(isempty(k))
        str = str(1:end-3);
%     end
end