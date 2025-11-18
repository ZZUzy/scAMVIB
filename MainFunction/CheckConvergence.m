function   [T,done] =CheckConvergence (T,prm,LoopCounter,change)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

% fprintf('Change is %d\n',change);
% fprintf('LoopCounter is %d\n\n',LoopCounter);

done =0;

if change <= prm.SmallChange
    %         T{w}.convergeStr = sprintf('sIB converged in iteration %d with change=%d',LoopCounter(w),change(w));
    done= 1;
elseif LoopCounter>= prm.LoopLimit
    %         T{w}.convergeStr = sprintf('sIB did NOT converge (change=%d), stop due to LoopLimit (LoopCounter=%d)',change(w),LoopCounter(w));
    done = 1;
end

end

