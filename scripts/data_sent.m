close all;
clc;
sent=[];
for i=1:length(Date)
    for j=1:length(Sent)
        sent(cast(Sent(j)*5,'uint8'),j) = Close(j);
    end
end

contourf(sent)
