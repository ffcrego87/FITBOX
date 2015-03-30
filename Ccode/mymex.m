%UNCOMMENT WHERE RELEVANT AND REPLACE DIRECTORIES

sprintf('subset...\n')
% %mex -LC:\ILOG\CPLEX112\lib\x86_windows_vs2008\stat_mda -lcplex112 -IC:\ILOG\cplex112\include\ilcplex subsetmex.c
%mex -L/home/pcasau/cplex121/lib/x86-64_debian4.0_4.1/static_pic -lcplex -I/home/pcasau/cplex121/include/ilcplex subsetmex.c

sprintf('lps...\n')
%mex -I'C:\ILOG\cplex112\include\ilcplex' C:\ILOG\cplex112\examples\src\check.c 'C:\ILOG\CPLEX112\lib\x86_windows_vs2008\stat_mda\cplex112.lib' lpsmex.c
%mex -L/home/pcasau/cplex121/lib/x86-64_debian4.0_4.1/static_pic/ -lcplex -I/home/pcasau/cplex121/include/ilcplex lpsmex.c

%mex -LC:\ILOG\CPLEX112\lib\x86_windows_vs2008\stat_mda -lcplex112 -IC:\ILOG\cplex112\include\ilcplex lpsmex.c
sprintf('reduce...\n')
%mex -LC:\ILOG\CPLEX112\lib\x86_windows_vs2008\stat_mda -lcplex112 -IC:\ILOG\cplex112\include\ilcplex reducemex.c
%mex -g -IC:\ILOG\cplex75\include\ilcplex bpsmex.c C:\ILOG\cplex75\lib\msvc6\stat_md\cplex75.lib
%mex -L/home/pcasau/cplex121/lib/x86-64_debian4.0_4.1/static_pic/ -lcplex -I/home/pcasau/cplex121/include/ilcplex reducemex.c
