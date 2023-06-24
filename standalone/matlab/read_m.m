function a=read_m(name)
fname=[name,'.txt'];
h=fopen(fname,'r');
if h<0,
    error(['cannot open file ',fname])
end
m=next(1,1);
n=next(1,1);
a=next(m,n);
h=fclose(h);
if h<0,
    error(['cannot close file ',fname])
end

    function a=next(mm,nn)  
        [a,count]=fscanf(h,'%g',[mm,nn]);
        if count ~= mm*nn,
            error(['premature end of file ',fname])
        end % if
    end % next

end % read_m
