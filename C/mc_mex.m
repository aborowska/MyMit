% [s,msg,msgid] = mkdir('c:\work');
% if (isempty(msgid))
%     mkdir('c:\work')
% end
% cd c:\work
% % Copy the source code, yprime.c.
% % 
copyfile(fullfile(matlabroot,'extern','examples','mex','yprime.c'),'.','f');
% Build the MEX-file.
mex yprime.c

T=1;
Y=1:4;
yprime(T,Y)

%%%%%%%%%%
edit([matlabroot '/extern/examples/mex/arrayProduct.c']);
addpath('include/');
addpath('C/');

mex C/arrayProduct.c

s = 5; 
A = [1.5, 2, 9];
B = arrayProduct(s,A)
% It is good practice to validate the type of a MATLAB variable before
% calling a MEX-file. To test the input variable, inputArg, and convert it
% to double, if necessary, use this code.
inputArg = int16(A);
if ~strcmp(class(inputArg),'double')
    inputArg = double(inputArg);
end
B = arrayProduct(s,inputArg)


edit([matlabroot '/extern/examples/mex/mexfunction.c']);

mex C/mexfunction.c

[d,e,f,dupa] = mexfunction(s,A,inputArg,theta,Y,msg)


mex C/play.c

theta2 = floor(10*rand(1,10));
[theta_R1,theta_R2] = play(theta,theta2)

mex C/try_mh.c