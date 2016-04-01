[s,msg,msgid] = mkdir('\\uva.nl\dfs\ti-home\ovu09568\Desktop\work');
if (isempty(msgid))
    mkdir('\\uva.nl\dfs\ti-home\ovu09568\Desktop\work')
end
cd '\\uva.nl\dfs\ti-home\ovu09568\Desktop\work'


copyfile(fullfile(matlabroot,'extern','examples','mex',...
    'yprime.c'),'.','f');

mex yprime.c

 mex.getCompilerConfigurations