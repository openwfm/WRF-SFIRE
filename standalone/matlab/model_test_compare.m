steps=6
out=read_model_test_out('model_test_out.txt',steps);
if ~exist('out_ref','var'),
    out_ref=read_model_test_out('model_test_out.txt.ref',steps);
end
for i=1:steps,max(max(abs((out.d(i).lfn-out_ref.d(i).lfn)))),end