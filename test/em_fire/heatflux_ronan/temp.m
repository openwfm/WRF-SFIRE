function temp
t=readout('T');
times=char(readout('Times')');
f=readout('GRNHFX');
for i=1:size(t,4)
    time=times(i,:);
    figure(1)
    mesh(t(:,:,1,i))
    title(['T level 1 ',time],'Interpreter', 'none')
    print('-dpng',['T1_',time,'.png'])
    figure(2)
    mesh(f(:,:,i))
    title(['GRNHFX ',time],'Interpreter', 'none')
    print('-dpng',['GRNHFX_',time,'.png'])
    drawnow,pause(0.5)
end        
end