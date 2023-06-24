function a=read_model_test_out(file,n)
f=fopen(file,'r');
a.dx=next_mat(f);
a.dy=next_mat(f);
t=next_mat(f)
i=0;
while ~isempty(t),
    i=i+1
    if i>n,break,end
    a.d(i).t=t;
    a.d(i).lfn=next_mat(f);
    a.d(i).tign=next_mat(f);
    a.d(i).vx=next_mat(f);
    a.d(i).vy=next_mat(f);
    a.d(i).flux=next_mat(f);
    t=next_mat(f)
end
fclose(f);
