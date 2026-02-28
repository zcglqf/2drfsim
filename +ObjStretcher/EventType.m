classdef EventType < uint8
    enumeration
        RF          (0)
        GX          (1)
        GY          (2)
        GZ          (3)
        ACQ         (4)
        EVENTMAX    (5)
    end

    methods
        function tf = lt(a,b)
            if isa(a,'ObjStretcher.EventType') && ...
               isa(b,'ObjStretcher.EventType')
                tf = uint8(a) < uint8(b);
            else
                error('ObjStretcher:EventType:TypeMismatch', ...
                      'Comparison must be between EventType values.');
            end
        end
    end
end