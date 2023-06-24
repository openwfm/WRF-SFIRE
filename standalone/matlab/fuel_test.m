function fuel_frac=fuel_test
for i=1:4,figure(i),end
%input('Posdslkjflkjdfition figure windows and press Enter >');
fd=[1,1];
[c1,c2]=ndgrid([0,fd(1)],[0,fd(2)]); % corners of the mesh cell
maxerr=0;

%tnow=1;
%fuel_time=1;
display('test1')
lfn0=[-1    -1; 1    1] % to debug fortran
tnow=2
%lfn0=   [0.4192792     -1.9766893E-02;   6.310914       5.983462 ]
fuel_time= 8.235294 ;
f_debug=[];
n=2^(9+1);

    tign = [1.000000      1.000000    ;
           2.000000       2.000000    ]
    f1=fuel_burnt(lfn0,tign,tnow,fd,fuel_time) 
    display('fuel_burnt result')
    f2=fuel_burnt_fd(lfn0,tign,tnow,fd,fuel_time) 
    display('fuel_burnt_fd result')
    fq=fuel_burnt_quad(lfn0,tign,tnow,fd,fuel_time,n)
    display('fuel_burnt_quad result')
    test1_ls=f1-fq       
    display('fuel_burnt error')
    test1_fd=f2-fq       
    display('fuel_burnt error')

    display('  ')
    
    display('test2')


lfn0=[-2    -2; -1    -1] % to debug fortran
tnow=3
%lfn0=   [0.4192792     -1.9766893E-02;   6.310914       5.983462 ]
fuel_time= 8.235294 ;
f_debug=[];
n=2^(9+1);

    tign = [1.000000      1.000000    ;
           2.000000       2.000000    ]
    f1=fuel_burnt(lfn0,tign,tnow,fd,fuel_time) 
    display('fuel_burnt result')
    f2=fuel_burnt_fd(lfn0,tign,tnow,fd,fuel_time) 
    display('fuel_burnt_fd result')
    fq=fuel_burnt_quad(lfn0,tign,tnow,fd,fuel_time,n)
    display('fuel_burnt_quad result')
    test2_ls=f1-fq       
    display('fuel_burnt error')
    test2_fd=f2-fq       
    display('fuel_burnt error')

    
display('test 3')

lfn0=[-2    -2; -4    -4] % to debug fortran
tnow=10
%lfn0=   [0.4192792     -1.9766893E-02;   6.310914       5.983462 ]
fuel_time= 8.235294 ;
f_debug=[];
n=2^(9+1);

    tign = [6.000000      6.000000    ;
           4.000000       4.000000    ]
    f1=fuel_burnt(lfn0,tign,tnow,fd,fuel_time) 
    display('fuel_burnt result')
    f2=fuel_burnt_fd(lfn0,tign,tnow,fd,fuel_time) 
    display('fuel_burnt_fd result')
    fq=fuel_burnt_quad(lfn0,tign,tnow,fd,fuel_time,n)
    display('fuel_burnt_quad result')
    test3_ls=f1-fq       
    display('fuel_burnt error')
    test3_fd=f2-fq       
    display('fuel_burnt error')

    
    
    
lfn0=[-3    -1; -1    1] % to debug fortran
tnow=10
%lfn0=   [0.4192792     -1.9766893E-02;   6.310914       5.983462 ]
fuel_time= 8.235294 ;
f_debug=[];
n=2^(9+1);

    tign = [4.000000      8.000000    ;
           8.000000       10.000000    ]
    f1=fuel_burnt(lfn0,tign,tnow,fd,fuel_time) 
    display('fuel_burnt result')
    f2=fuel_burnt_fd(lfn0,tign,tnow,fd,fuel_time) 
    display('fuel_burnt_fd result')
    fq=fuel_burnt_quad(lfn0,tign,tnow,fd,fuel_time,n)
    display('fuel_burnt_quad result')
    test4_ls=f1-fq       
    display('fuel_burnt error')
    test4_fd=f2-fq       
    display('fuel_burnt error')

lfn0=[-1    1; 1    3] % to debug fortran
tnow=2
%lfn0=   [0.4192792     -1.9766893E-02;   6.310914       5.983462 ]
fuel_time= 8.235294 ;
f_debug=[];
n=2^(9+1);

    tign = [1.000000      2.000000    ;
           2.000000       2.000000    ]
    f1=fuel_burnt(lfn0,tign,tnow,fd,fuel_time) 
    display('fuel_burnt result')
    f2=fuel_burnt_fd(lfn0,tign,tnow,fd,fuel_time) 
    display('fuel_burnt_fd result')
    fq=fuel_burnt_quad(lfn0,tign,tnow,fd,fuel_time,n)
    display('fuel_burnt_quad result')
    test5_ls=f1-fq       
    display('fuel_burnt error')
    test5_fd=f2-fq       
    display('fuel_burnt error')
    test1_ls
    test1_fd
   
    test2_ls
    test2_fd
    
    test3_ls
    test3_fd
    
    test4_ls
    test4_fd
    
    test5_ls
    test5_fd
   
    
    
    
end


function fuel_frac=fuel_burnt_quad(lfn,tign,tnow,fd,fuel_time,n)
% compute the same numerically for comparison
% if fireline passes through the cell then
% lfn and tign must be linear
%otherwise the result will differ
% evaluate by quadrature on n by n points
% the mesh
g=(0.5+[0:n-1])/n;
[q1,q2]=meshgrid(fd(1)*g,fd(2)*g); % quad nodes
[c1,c2]=meshgrid([0,fd(1)],[0,fd(2)]); % corners of the mesh cell
% T and L on quad nodes
T = interp2(c1,c2,tign,q1,q2)-tnow;
L = interp2(c1,c2,lfn,q1,q2);
% the integrand
f = (1 - exp((L<0).*T./fuel_time)).*(L<0);
% integrate
fuel_frac = sum(f(:))/(n*n);
end

