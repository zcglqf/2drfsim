function samples = metaToSamples(metaArray)
%ObjStretcher.metaToSamples Convert EventMetaData[] to EventSampleData[].
%
% Each EventMetaData fragment generates two time points:
%   (StartTime, StartOrEnd=true , Type, AverageValue)
%   (EndTime  , StartOrEnd=false, Type, AverageValue)

n = numel(metaArray);
samples(2*n,1) = ObjStretcher.EventSampleData(); % preallocate

for i = 1:n
    m = metaArray(i);

    samples(2*i-1) = ObjStretcher.EventSampleData( ...
        't', m.StartTime, ...
        'startOrEnd', true, ...
        'type', m.Type, ...
        'v', m.AverageValue);

    samples(2*i) = ObjStretcher.EventSampleData( ...
        't', m.EndTime, ...
        'startOrEnd', false, ...
        'type', m.Type, ...
        'v', m.AverageValue);
end

end