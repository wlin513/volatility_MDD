function [ out ] = arctg( innum )

    if innum < 0 | innum > 1
        error('inverse logit only defined for numbers between 0 and 1');
    end
   % out=1/pi*(atan(innum)+pi/2);
out=atan(pi*(2*innum-1)/2);
end

