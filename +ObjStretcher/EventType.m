classdef EventType < uint8
    enumeration
        RF          (0)
        GRADIENTM   (3)
        GRADIENTP   (4)
        GRADIENTS   (5)
        ACQ         (6)
        EVENTMAX    (7)
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