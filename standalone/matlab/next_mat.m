   function mat=next_mat(f)
    %mat=reshape(a(k:k+mm*nn),mm,nn);
    [sz,count]=fscanf(f,'%g',[1 2]);
    switch count
        case 0,
            mat=[];
        case 1
            error('not enough terms in size')
        otherwise
        %fprintf('reading matrix size %g %g\n',sz)
        if any(any(sz~=round(sz) | sz < 0)),
            error('bad size')
        end
        [mat,count]=fscanf(f,'%g',sz);
        if(count<sz(1)*sz(2)),
            error('not enough terms in matrix')
        end
    end
    end
