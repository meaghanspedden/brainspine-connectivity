
function sensornumberout = get_sensornumber(sensname)
%input formatted as
% 'R8'

if contains(sensname, 'NC')
    sensnums=1:8;
elseif contains(sensname,'W')
    sensnums=9:16;
elseif contains(sensname, 'LB')
    sensnums=17:24;
elseif contains(sensname, 'R')
    sensnums=25:32;
elseif contains(sensname, 'BR')
    sensnums=33:40;
elseif contains(sensname, 'SW')
    sensnums=41:48;
elseif contains(sensname, 'DB')
    sensnums=49:56;
elseif contains(sensname, 'G')
    sensnums=57:64;
end

inputnum=str2double(sensname(end));
sensornumber=int2str(sensnums(inputnum));

sensornumberout = {['X' sensornumber], ['Y' sensornumber], ['Z' sensornumber]};



