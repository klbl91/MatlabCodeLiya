function [A,Eo] = BackAllFWD(A,varargin)

j = 0;
for i=1:length(A)
    fprintf('Doing %d/%d\n',i,length(A));
    
    for iDrop = 1:length(A(i).Drops)        
        try
            temp = KalmanBack(A(i),'iDrop',iDrop,varargin{:});
            
            j = j+1;
            B(j).Result = temp;
        catch
            fprintf('Error occured: %s\n',lasterr);
        end
    end
    
end

A = B;


